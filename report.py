
# FASE RESULTS PLOTTING AND REPORTING
# TEMPORARILY IMPLEMENTED AS A STANDALONE SCRIPT

from typing import List, Any
import os
import configparser
from functools import reduce

from src.analysis.experiment import Sample
from src.analysis.counts import run_feature_counts, get_gene_counts_from_tsv, get_tmm_cpm_from_gene_counts
from src.reporting.process_results import load_fase_results
from src.reporting.plot import plot_transcript, plot_splice_rate
from src.reporting.generate_html_report import generate_html_report


REFERENCE_GENOME_GTF_PATH = ("/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCm39/"
                             "gencode_M27_primary/gencode.vM27.primary_assembly.annotation.gtf")

BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/hons_PD1KO/sorted"
BAM_SUFFIX = "_Aligned_Sorted.out.bam"

FASE_RESULTS_PATH = "/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4 SplDC2020/transmembrane region.csv"

FEATURECOUNTS_EXECUTABLE = "/Users/aasho2/opt/anaconda3/envs/bbmap/bin/featureCounts"
PRIMARY_ALIGNMENT_ONLY = True
N_THREADS = 24

RUN_CONFIG_PATH = "/Users/aasho2/Projects/FASE_test/fase_run.config"
OUTPUT_DIR = "OUTPUT/Report"


def _check_inputs_are_valid(
    bam_file_absolute_paths: List[str],
    bam_files_dir: str,
    bam_suffix: str,
    fase_results_path: str,
    output_dir: str,
) -> None:

    # Check BAM files
    if len(bam_file_absolute_paths) == 0:
        raise ValueError(f"No BAM files detected in {bam_files_dir} with suffix {bam_suffix}")

    # Check FASE results file
    if not os.path.exists(fase_results_path):
        raise FileNotFoundError("The specified FASE results file was not found")

    # Check output directory (but do not create it yet)
    if not os.path.isdir(output_dir) and os.path.exists(output_dir):
        raise IsADirectoryError("The specified output directory is a file path")


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def main(
    bam_files_dir: str,
    bam_suffix: str,
    fase_results_path: str,
    output_dir: str
) -> None:

    # ---------------------------------------------------------------------------------------------------------------- #

    # Get BAM file paths
    _suffix_length = len(bam_suffix)
    bam_file_absolute_paths = [os.path.join(
        bam_files_dir, p
    ) for p in os.listdir(bam_files_dir) if p[-_suffix_length:] == bam_suffix]

    # Ensure output_dir is an absolute path
    output_dir_absolute = os.path.abspath(output_dir)

    # Proceed only if inputs pass validation
    _check_inputs_are_valid(bam_file_absolute_paths, bam_files_dir, bam_suffix, fase_results_path, output_dir_absolute)

    # Create output directory if required
    if not os.path.isdir(output_dir_absolute):
        os.makedirs(output_dir_absolute)

    # Load samples
    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, BAM_SUFFIX)
        samples[sample.name] = sample

    # ---------------------------------------------------------------------------------------------------------------- #

    # Load run config
    run_config = configparser.ConfigParser()
    run_config.read(RUN_CONFIG_PATH)

    # TODO: Refactor pointless assertions - missing values will trigger KeyError before reaching assert statement
    run_name = run_config["GENERAL"]["name"]
    assert run_name, f"A run name is required (attribute 'name' in section [GENERAL])"

    feature_name = run_config["FEATURE"]["match"]
    assert run_name, f"A feature pattern match is required (attribute 'match' in section [FEATURE])"

    sample_groups = {key: val.split(" ") for key, val in run_config["SAMPLES"].items() if val}

    assert all(
        [sample_name in samples for sample_name in flatten_nested_lists(list(sample_groups.values()))]
    ), "Sample names in run config do not match sample names parsed from BAM files"
    assert all(
        [sample_name in flatten_nested_lists(list(sample_groups.values())) for sample_name in samples]
    ), "Sample names parsed from BAM files are not all present in the SAMPLES section of run config"

    group_name_by_sample = {}
    for group_name, sample_names in sample_groups.items():
        for sample_name in sample_names:
            group_name_by_sample[sample_name] = group_name

    color_by_group_name = dict(run_config["COLORS"])

    assert all(
        [group_name in sample_groups.keys() for group_name in color_by_group_name.keys()]
    ), "A group name is defined in the COLORS section of run config that is not defined in the SAMPLES section"
    assert all(
        [group_name in color_by_group_name.keys() for group_name in sample_groups.keys()]
    ), "A group name is defined in the SAMPLES section of run config that is not defined in the COLORS section"

    # ---------------------------------------------------------------------------------------------------------------- #

    # Get TMM CPM
    gene_counts_path = os.path.join(output_dir_absolute, "gene_counts.tsv")

    if not os.path.exists(gene_counts_path):
        print("Running feature counts...")
        run_feature_counts(
            FEATURECOUNTS_EXECUTABLE,
            BAM_FILES_DIR,
            BAM_SUFFIX,
            REFERENCE_GENOME_GTF_PATH,
            gene_counts_path,
            paired_end_reads=PRIMARY_ALIGNMENT_ONLY,
            threads=N_THREADS
        )
        print("... done")

    print("Reading gene counts...")

    raw_gene_counts = get_gene_counts_from_tsv(
        gene_counts_path,
        bam_suffix=BAM_SUFFIX
    )

    print("... done")

    print("Applying TMM normalization to gene counts...")

    norm_gene_counts = get_tmm_cpm_from_gene_counts(raw_gene_counts)

    print("... done")

    # ---------------------------------------------------------------------------------------------------------------- #

    # Load FASE results
    fase_results = load_fase_results(
        fase_results_path,
        samples,
        min_total_number_occurrences_across_all_samples=999999,
        #     force_include_gene_names = ["Clec9a"]
        #     force_include_gene_names = ["H2-Q7"]
        #     force_include_gene_names = ["Cd86"]
        #     force_include_gene_names = ["Treml4"]
        #     force_include_gene_names = ["Ifnar2", "Ifnlr1"]
        force_include_gene_names=["Tnfrsf9", "Cd69"]
        #     force_include_gene_names = ["Pdcd1", "Cd274", "Pdcd1lg2"]
    )

    print("Generating figures for report...")

    plots = {
        fase_result.transcript_id:
        {
            "expression": plot_splice_rate(
                fase_result,
                norm_gene_counts,
                sample_groups,
                group_name_by_sample,
                color_by_group_name
            ),
            "splice": plot_transcript(
                fase_result,
                draw_junctions_with_min_n_occurrences=2,
                show_main_title=False
            )
        }
        for fase_result in fase_results
    }

    print("... finished")

    print("Generating report...")

    report_html = generate_html_report(
        run_name,
        feature_name,
        fase_results,
        plots
    )

    with open(os.path.join(output_dir_absolute, "report.html"), "w") as output_file:
        output_file.write(report_html)

    print("...finished")


if __name__ == "__main__":

    main(
        BAM_FILES_DIR,
        BAM_SUFFIX,
        FASE_RESULTS_PATH,
        OUTPUT_DIR
    )
