
# FASE RESULTS PLOTTING AND REPORTING
# TEMPORARILY IMPLEMENTED AS A STANDALONE SCRIPT

from typing import List, Dict, Any
import os
from pandas import read_csv as pd_read_csv
from functools import reduce

from src.config.parse_config import FaseRunConfig
from src.analysis.core import name_output_file
from src.analysis.experiment import Sample
from src.analysis.counts import run_feature_counts, get_gene_counts_from_tsv, get_tmm_cpm_from_gene_counts
from src.reporting.process_results import load_fase_results
from src.reporting.plot import plot_transcript, plot_splice_rate
from src.reporting.generate_html_report import generate_html_report


GENE_COUNTS_DEFAULT_FILE_NAME = "gene_counts.tsv"
TMM_NORM_GENE_COUNTS_DEFAULT_FILE_NAME = "tmm_norm_gene_counts.tsv"


def _check_inputs_are_valid(
    bam_file_absolute_paths: List[str],
    bam_files_dir: str,
    bam_ending: str,
    fase_results_path: str,
    output_dir: str,
) -> None:

    # TODO: Refactor to just take run_config and check attributes (are these checks redundant if __init__ passed?)

    # Check BAM files
    if len(bam_file_absolute_paths) == 0:
        raise ValueError(f"No BAM files detected in {bam_files_dir} with ending {bam_ending}")

    # Check FASE results file
    if not os.path.exists(fase_results_path):
        raise FileNotFoundError(f"The specified FASE results file was not found: {fase_results_path}")

    # Check output directory (but do not create it yet)
    if not os.path.isdir(output_dir) and os.path.exists(output_dir):
        raise IsADirectoryError("The specified output directory is a file path")


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def create_report(
    run_config: FaseRunConfig
    # bam_files_dir: str,
    # bam_ending: str,
    # fase_results_path: str,
    # output_dir: str,
    # max_n_plotted: Union[None, int] = None
) -> None:

    # ---------------------------------------------------------------------------------------------------------------- #

    output_dir = os.path.join(run_config.output_path, "REPORT")
    _results_file_name_no_extension = name_output_file(run_config.run_name, run_config.feature_name)
    fase_results_path = os.path.join(run_config.output_path, f"{_results_file_name_no_extension}.csv")

    # Get BAM file paths
    _ending_length = len(run_config.bam_ending)
    bam_file_absolute_paths = [os.path.join(
        run_config.input_path, p
    ) for p in os.listdir(run_config.input_path) if p[-_ending_length:] == run_config.bam_ending]

    # Ensure output_dir is an absolute path
    output_dir_absolute = os.path.abspath(output_dir)

    # Proceed only if inputs pass validation
    _check_inputs_are_valid(
        bam_file_absolute_paths,
        run_config.input_path,
        run_config.bam_ending,
        fase_results_path,
        output_dir_absolute
    )

    # Create output directory if required
    if not os.path.isdir(output_dir_absolute):
        os.makedirs(output_dir_absolute)

    # Load samples
    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, run_config.bam_ending)
        samples[sample.name] = sample

    # Get TMM CPM
    gene_counts_path = run_config.report_gene_count_matrix \
        if run_config.report_gene_count_matrix is not None \
        else os.path.join(output_dir_absolute, GENE_COUNTS_DEFAULT_FILE_NAME)

    if not os.path.exists(gene_counts_path):
        print("Running feature counts...")
        run_feature_counts(
            run_config.report_feature_counts_executable,
            run_config.input_path,
            run_config.bam_ending,
            run_config.genome,
            gene_counts_path,
            paired_end_reads=run_config.paired_end_reads,
            primary_alignment_only=run_config.primary_alignment_only,
            threads=run_config.report_n_threads
        )
        print("... done")
    else:
        print("An existing gene count matrix has been supplied.\n")

    tmm_norm_gene_counts_path = os.path.join(output_dir_absolute, TMM_NORM_GENE_COUNTS_DEFAULT_FILE_NAME)

    # Load TMM norm counts if saved, otherwise read raw counts and perform normalization
    if os.path.isfile(tmm_norm_gene_counts_path):

        print("Reading TMM normalized gene counts...")

        norm_gene_counts = pd_read_csv(
            tmm_norm_gene_counts_path,
            sep="\t"
        )
        norm_gene_counts.set_index("Geneid", inplace=True)

        print("...done")

    else:

        print("Reading gene counts...")

        raw_gene_counts = get_gene_counts_from_tsv(
            gene_counts_path,
            bam_ending=run_config.bam_ending
        )

        print("... done")

        print("Applying TMM normalization to gene counts...")

        norm_gene_counts = get_tmm_cpm_from_gene_counts(
            raw_gene_counts,
            threads=run_config.n_threads
        )

        # Write TMM norm CPM to file if specified to do so in config
        if run_config.report_save_normalized_gene_counts:
            norm_gene_counts["Geneid"] = norm_gene_counts.index
            norm_gene_counts.to_csv(
                tmm_norm_gene_counts_path,
                sep="\t",
                index=False
            )

        print("... done")

    # ---------------------------------------------------------------------------------------------------------------- #

    # Load FASE results
    fase_results = load_fase_results(
        fase_results_path,
        samples,
        rank_by=run_config.rank_results_by
        # min_total_number_occurrences_across_all_samples=999999,
        #     force_include_gene_names = ["Clec9a"]
        #     force_include_gene_names = ["H2-Q7"]
        #     force_include_gene_names = ["Cd86"]
        #     force_include_gene_names = ["Treml4"]
        #     force_include_gene_names = ["Ifnar2", "Ifnlr1"]
        # force_include_gene_names=["Tnfrsf9", "Cd69"]
        #     force_include_gene_names = ["Pdcd1", "Cd274", "Pdcd1lg2"]
    )

    print("Generating figures for report...")

    plots: Dict[str, Dict[str, str]] = {}
    _current = 1
    _total = len(fase_results) if run_config.report_max_n_plotted is None \
        else min(len(fase_results), run_config.report_max_n_plotted)
    print()  # Print blank line to be consumed by first one-line-up ansi code
    for fase_result in fase_results:

        if run_config.report_max_n_plotted is not None and _current > run_config.report_max_n_plotted:
            break

        # TODO: Fix so that entire line is overwritten
        _temp_solution_right_spacing = "          "
        print(
            "\x1b[A" + f"Progress: [{_current}/{_total}] {fase_result.gene_name} / {fase_result.transcript_id} " +
            f"feature {fase_result.feature_number} of {fase_result.total_features_in_transcript}" +
            _temp_solution_right_spacing
        )
        _current += 1

        plots[f"{fase_result.transcript_id}-{fase_result.feature_number}"] = {
            "expression": plot_splice_rate(
                fase_result,
                norm_gene_counts,
                run_config.sample_groups,
                run_config.group_name_by_sample_name,
                run_config.color_by_group_name,
                show_main_title=False
            ),
            "splice": plot_transcript(
                fase_result,
                # TODO: Set as parameter in run config REPORT section
                draw_junctions_with_min_n_occurrences=run_config.report_draw_junctions_min_count,
                show_main_title=False
            )
        }

    print("... finished")

    print("Generating report...")

    report_html = generate_html_report(
        run_config.run_name,
        run_config.feature_name,
        fase_results,
        plots
    )

    output_path = os.path.join(
        output_dir_absolute, f"Report - {name_output_file(run_config.report_name, run_config.feature_name)}.html"
    )

    with open(output_path, "w") as output_file:
        output_file.write(report_html)

    print(f"Report written to {output_path}")
    print("...finished")
