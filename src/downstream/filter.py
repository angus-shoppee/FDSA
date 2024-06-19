# Feature-Directed Splice Analysis
# Copyright (C) 2024 Angus Shoppee
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from typing import List, Tuple
import os
import subprocess

from config.parse_config import ProgramRunConfig
from analysis.core import name_output, Sample
from reporting.process_results import load_fdsa_results


SAM_FLAG_REVERSE_STRAND = 16  # If this flag is present in a read, it is mapped to the reverse strand

FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME = "FILTERED_BAM"


def get_gene_start_and_end(exon_positions: List[Tuple[int, int]]) -> Tuple[int, int]:

    flattened_coordinates = [coordinate for exon_position in exon_positions for coordinate in exon_position]

    return min(flattened_coordinates), max(flattened_coordinates)


def generate_filtered_bam_files(
    run_config: ProgramRunConfig,
    check_exit_code: bool = True
) -> None:

    output_dir = os.path.abspath(
        os.path.join(
            run_config.output_path,
            FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME,
            f"{name_output(run_config.run_name, run_config.feature_name)}"
        )
    )
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_files_dir = os.path.abspath(run_config.input_path)

    samtools = run_config.filter_samtools_executable
    if not os.path.isfile(samtools):
        raise ValueError(f"The supplied samtools executable path is not a valid file")

    fdsa_results_path = os.path.join(
        run_config.output_path,
        f"{name_output(run_config.run_name, run_config.feature_name)}.csv"
    )

    # Get BAM file paths
    _ending_length = len(run_config.bam_ending)
    bam_file_absolute_paths = [os.path.join(
        bam_files_dir, p
    ) for p in os.listdir(run_config.input_path) if p[-_ending_length:] == run_config.bam_ending]

    # Define samples
    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, run_config.bam_ending)
        samples[sample.name] = sample

    fdsa_results = load_fdsa_results(
        fdsa_results_path,
        samples,
        min_total_number_occurrences_across_all_samples=run_config.filter_min_total_n_occurrences_across_all_samples,
        min_per_sample_occurrences_number_occurrences=run_config.filter_min_n_occurrences_in_sample,
        min_per_sample_occurrences_in_at_least_n_samples=run_config.filter_occurrences_in_at_least_n_samples
    )

    # REMOVED: Running samtools view once with a BED file specifying multiple regions is far slower (why?) than
    #          running individually for each region, despite producing the same output

    # # Create BED file specifying regions
    #
    # bed_file_path = os.path.join(output_dir, "regions.bed")
    #
    # with open(bed_file_path, "w") as bed_file:
    #
    #     for result in fdsa_results:
    #
    #         span = get_gene_start_and_end(result.exon_positions)
    #         # genomic_region = f"{result.chromosome}:{span[0]}-{span[1]}"
    #
    #         bed_file.write(f"{result.chromosome}\t{span[0]}\t{span[1]}\n")

    print("Generating filtered BAM files...")

    _right_padding = "            "
    _n_results = len(fdsa_results)

    for sample in samples.values():

        print(f"Extracting reads from: {sample.name}{sample.bam_ending}")
        print()  # Include blank line to be consumed by first one-line-up ansi code

        unsorted_sam_out_path = os.path.join(
            output_dir,
            f"{sample.name}{sample.bam_ending.replace(run_config.bam_ending, '.temp.sam')}"
        )

        with open(unsorted_sam_out_path, "w") as out_file:

            # REMOVED: BED file method

            # split_cmd = [
            #     samtools, "view",
            #     # "-b",
            #     "--with-header",
            #     "--target-file", bed_file_path,
            #     sample.bam_path
            # ]
            #
            # print("[DEBUG] Trying:")
            # print(" ".join(split_cmd))
            #
            # _process = subprocess.run(
            #     split_cmd,
            #     cwd=output_dir,
            #     encoding="utf8",
            #     check=check_exit_code,
            #     stdout=out_file
            # )

            # END BED file method

            # First, extract the header

            split_cmd = [
                samtools, "view", sample.bam_path, "--header-only"
            ]

            _process = subprocess.run(
                split_cmd,
                cwd=output_dir,
                encoding="utf8",
                check=check_exit_code,
                stdout=out_file
            )

            # Then, write reads to output SAM file for each qualifying gene

            for i, result in enumerate(fdsa_results):

                print("\x1b[A" + f"[{i + 1} / {_n_results}] {result.gene_name}" + _right_padding)

                span = get_gene_start_and_end(result.exon_positions)
                genomic_region = f"{result.chromosome}:{span[0]}-{span[1]}"

                if result.strand == "+":
                    strand_specifier = ["--exclude-flags", str(SAM_FLAG_REVERSE_STRAND)]
                elif result.strand == "-":
                    strand_specifier = ["--require-flags", str(SAM_FLAG_REVERSE_STRAND)]
                else:
                    raise ValueError(f"Encountered FDSA output for gene {result.gene_name} with invalid strand "
                                     f"(\"{result.strand}\"). Expected either \"+\" or \"-\"")

                # split_cmd = [
                #     samtools, "view", sample.bam_path, genomic_region
                # ] + strand_specifier
                split_cmd = [
                    samtools, "view", sample.bam_path, genomic_region
                ]
                # TODO: Strand selection has been removed for now - unsure whether mate reads were handled correctly

                if run_config.filter_unique_mapping_only:
                    split_cmd += ["--min-MQ", str(run_config.mapq_for_unique_mapping)]

                _process = subprocess.run(
                    split_cmd,
                    cwd=output_dir,
                    encoding="utf8",
                    check=check_exit_code,
                    stdout=subprocess.PIPE
                )

                if _process.stdout:

                    out_file.write(_process.stdout)

            # Allow file to close before proceeding to final step

        # Sort reads and output as BAM file

        sorted_out_path = unsorted_sam_out_path.replace(".temp.sam", run_config.bam_ending)
        # sorted_out_path = unsorted_sam_out_path.replace(".temp.sam", run_config.bam_ending).replace(".bam", ".sam")

        print("Sorting...")

        split_cmd = [
            samtools, "sort",
            "--output-fmt", "BAM",
            "--threads", str(run_config.n_threads),
            "-o", sorted_out_path,
            unsorted_sam_out_path
        ]

        _process = subprocess.run(
            split_cmd,
            cwd=output_dir,
            encoding="utf8",
            check=check_exit_code,
            stderr=subprocess.DEVNULL
        )

        print("... done")

        # Index the final BAM file

        print("Indexing...")

        split_cmd = [
            samtools, "index",
            sorted_out_path
        ]

        _process = subprocess.run(
            split_cmd,
            cwd=output_dir,
            encoding="utf8",
            check=check_exit_code
        )

        print("... done")

        # Delete unsorted temp file
        os.remove(unsorted_sam_out_path)

    print("... finished\n")
