
from typing import List, Tuple
import os
import subprocess

from src.config.parse_config import FaseInternalConfig, FaseRunConfig
from src.analysis.core import name_output_file, Sample
from src.reporting.process_results import load_fase_results


SAM_FLAG_REVERSE_STRAND = 16  # If this flag is present in a read, it is mapped to the reverse strand

FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME = "FILTERED_BAM"


def get_gene_start_and_end(exon_positions: List[Tuple[int, int]]) -> Tuple[int, int]:

    flattened_coordinates = [coordinate for exon_position in exon_positions for coordinate in exon_position]

    return min(flattened_coordinates), max(flattened_coordinates)


def generate_filtered_bam_files(
    run_config: FaseRunConfig,
    check_exit_code: bool = True
) -> None:

    output_dir = os.path.abspath(
        os.path.join(
            run_config.output_path,
            FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME,
            f"{name_output_file(run_config.run_name, run_config.feature_name)}"
        )
    )
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files_dir = os.path.abspath(run_config.input_path)

    samtools = run_config.filter_samtools_executable
    if not os.path.isfile(samtools):
        raise ValueError(f"The supplied samtools executable path is not a valid file")

    # def load_fase_results(
    #     fase_results_path: str,
    #     samples: Dict[str, Sample],
    #     rank_by: str = "frequency",
    #     force_gene_names: Union[None, List[str]] = None,
    #     min_total_number_occurrences_across_all_samples: int = 1,
    #     min_per_sample_occurrences_number_occurrences: int = 0,
    #     min_per_sample_occurrences_in_at_least_n_samples: int = 0,
    #     fase_result_frequency_column_prefix: str = FASE_RESULT_FREQUENCY_COLUMN_PREFIX
    # ) -> List[FaseResult]:

    fase_results_path = os.path.join(
        run_config.output_path,
        f"{name_output_file(run_config.run_name, run_config.feature_name)}.csv"
    )

    # Get BAM file paths
    _ending_length = len(run_config.bam_ending)
    bam_file_absolute_paths = [os.path.join(
        run_config.input_path, p
    ) for p in os.listdir(run_config.input_path) if p[-_ending_length:] == run_config.bam_ending]

    # Define samples
    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, run_config.bam_ending)
        samples[sample.name] = sample

    fase_results = load_fase_results(
        fase_results_path,
        samples,
        min_total_number_occurrences_across_all_samples=run_config.filter_min_total_n_occurrences_across_all_samples,
        min_per_sample_occurrences_number_occurrences=run_config.filter_min_n_occurrences_in_sample,
        min_per_sample_occurrences_in_at_least_n_samples=run_config.filter_occurrences_in_at_least_n_samples
    )

    print("Generating filtered BAM files...")

    _right_padding = "            "
    _n_results = len(fase_results)

    for sample in samples.values():

        print(f"Extracting reads from: {sample.name}{sample.bam_ending}")
        print()  # Include blank line to be consumed by first one-line-up ansi code

        sam_out_path = os.path.join(output_dir, f"{sample.name}{sample.bam_ending.replace('.bam', '.sam')}")

        mapq_specifier = f" --min-MQ {run_config.mapq_for_unique_mapping}" \
            if run_config.filter_unique_mapping_only else ""  # Leading space if not blank

        with open(sam_out_path, "w") as out_file:

            # First, extract the header

            full_cmd = f"{samtools} view {sample.bam_path} --header-only"

            process = subprocess.run(
                full_cmd.split(" "),
                cwd=output_dir,
                encoding="utf8",
                check=check_exit_code,
                stdout=out_file
            )

            # Then, write reads to output SAM file for each qualifying gene

            for i, result in enumerate(fase_results):
                print("\x1b[A" + f"[{i + 1} / {_n_results}] {result.gene_name}" + _right_padding)

                span = get_gene_start_and_end(result.exon_positions)
                genomic_region = f"{result.chromosome}:{span[0]}-{span[1]}"

                out_file.write("\n")  # Add a newline to separate appended content

                if result.strand == "+":
                    strand_specifier = f"--exclude-flags {SAM_FLAG_REVERSE_STRAND}"
                elif result.strand == "-":
                    strand_specifier = f"--require-flags {SAM_FLAG_REVERSE_STRAND}"
                else:
                    raise ValueError(f"Encountered FASE output for gene {result.gene_name} with invalid strand "
                                     f"(\"{result.strand}\"). Expected either \"+\" or \"-\"")

                # Lack of space between {strand_specifier} and {mapq_specifier} is intentional
                full_cmd = f"{samtools} view {sample.bam_path} {genomic_region} {strand_specifier}{mapq_specifier}"

                process = subprocess.run(
                    full_cmd.split(" "),
                    cwd=output_dir,
                    encoding="utf8",
                    check=check_exit_code,
                    stdout=out_file
                )

            # Allow file to close before proceeding to final step

        # Finally, sort output and return to BAM format

        # TODO: Enable "--output-fmt BAM" when finished testing

        # full_cmd = f"{samtools} sort --output-fmt BAM --threads {run_config.n_threads} {sam_out_path}"
        full_cmd = f"{samtools} sort --threads {run_config.n_threads} {sam_out_path}"

        process = subprocess.run(
            full_cmd.split(" "),
            cwd=output_dir,
            encoding="utf8",
            check=check_exit_code
        )

    print("... done\n")

    pass
