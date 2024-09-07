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


# FDSA RESULTS PLOTTING AND REPORTING

from typing import List, Dict, Any
import logging
import os
from pandas import read_csv as pd_read_csv
from functools import reduce

from config.parse_config import ProgramRunConfig
from analysis.core import name_output
from analysis.experiment import Sample
from analysis.counts import run_feature_counts, get_gene_counts_from_tsv, get_tmm_cpm_from_gene_counts
from reporting.process_results import load_fdsa_results
from reporting.plot import plot_transcript, plot_splice_rate
from reporting.generate_html_report import generate_html_report
from downstream.quant import STRINGTIE_OUTPUT_DIR_NAME, COMBINED_RESULTS_FILE_NAME


GENE_COUNTS_DEFAULT_FILE_NAME = "gene_counts.tsv"
TMM_NORM_GENE_COUNTS_DEFAULT_FILE_NAME = "tmm_norm_gene_counts.tsv"


logger = logging.getLogger(__name__)


def _check_inputs_are_valid(
    bam_file_absolute_paths: List[str],
    bam_files_dir: str,
    bam_ending: str,
    fdsa_results_path: str,
    output_dir: str,
) -> None:

    # TODO: Refactor to just take run_config and check attributes (are these checks redundant if __init__ passed?)

    # Check BAM files
    if len(bam_file_absolute_paths) == 0:
        raise ValueError(f"No BAM files detected in {bam_files_dir} with ending {bam_ending}")

    # Check FDSA results file
    if not os.path.exists(fdsa_results_path):
        raise FileNotFoundError(f"The specified FDSA results file was not found: {fdsa_results_path}")

    # Check output directory (but do not create it yet)
    if not os.path.isdir(output_dir) and os.path.exists(output_dir):
        raise IsADirectoryError("The specified output directory is a file path")


def flatten_nested_lists(nested_lists: List[List[Any]]) -> List[Any]:
    return reduce(lambda a, b: a + b, nested_lists)


def create_report(
    run_config: ProgramRunConfig
    # bam_files_dir: str,
    # bam_ending: str,
    # fdsa_results_path: str,
    # output_dir: str,
    # max_n_plotted: Union[None, int] = None
) -> None:

    # ---------------------------------------------------------------------------------------------------------------- #

    output_dir = os.path.join(run_config.output_path, "REPORT")
    _results_file_name_no_extension = name_output(run_config.run_name, run_config.feature_name)
    fdsa_results_path = os.path.join(run_config.output_path, f"{_results_file_name_no_extension}.csv")

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
        fdsa_results_path,
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
        logger.info("Running feature counts...")
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
        logger.info("... done")
    else:
        logger.info("An existing gene count matrix has been supplied.\n")

    tmm_norm_gene_counts_path = os.path.join(output_dir_absolute, TMM_NORM_GENE_COUNTS_DEFAULT_FILE_NAME)

    # Load TMM norm counts if saved, otherwise read raw counts and perform normalization
    if os.path.isfile(tmm_norm_gene_counts_path):

        logger.info("Reading TMM normalized gene counts...")

        norm_gene_counts = pd_read_csv(
            tmm_norm_gene_counts_path,
            sep="\t"
        )
        norm_gene_counts.set_index("Geneid", inplace=True)

        logger.info("...done")

    else:

        logger.info("Reading gene counts...")

        raw_gene_counts = get_gene_counts_from_tsv(
            gene_counts_path,
            bam_ending=run_config.bam_ending
        )

        logger.info("... done")

        logger.info("Applying TMM normalization to gene counts...")

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

        logger.info("... done")

    # ---------------------------------------------------------------------------------------------------------------- #

    # TODO: Revamp method for detecting quant mode output
    stringtie_results_path = os.path.join(
        os.path.abspath(run_config.output_path), STRINGTIE_OUTPUT_DIR_NAME, f"{COMBINED_RESULTS_FILE_NAME}.csv"
    )

    # Load FDSA results
    fdsa_results = load_fdsa_results(
        fdsa_results_path,
        samples,
        rank_by=run_config.rank_results_by,
        force_gene_names=None if run_config.report_genes is None else list(run_config.report_genes.keys()),
        min_total_number_occurrences_across_all_samples=run_config.report_min_total_n_occurrences_across_all_samples,
        min_per_sample_occurrences_number_occurrences=run_config.report_min_n_occurrences_in_sample,
        min_per_sample_occurrences_in_at_least_n_samples=run_config.report_occurrences_in_at_least_n_samples,
        stringtie_results_path=stringtie_results_path if os.path.isfile(stringtie_results_path) else None,
        stringtie_remove_sample_name_ending=run_config.bam_ending
    )

    logger.info("Generating figures for report...")

    plots: Dict[str, Dict[str, str]] = {}
    _current = 1
    _total = len(fdsa_results) if run_config.report_max_n_plotted is None \
        else min(len(fdsa_results), run_config.report_max_n_plotted)
    # logger.info()  # Print blank line to be consumed by first one-line-up ansi code
    for fdsa_result in fdsa_results:

        if run_config.report_max_n_plotted is not None and _current > run_config.report_max_n_plotted:
            break

        # _temp_solution_right_spacing = "            "
        logger.info(
            # "\x1b[A" +
            f"Progress: [{_current}/{_total}] {fdsa_result.gene_name} / {fdsa_result.transcript_id} " +
            f"feature {fdsa_result.feature_number} of {fdsa_result.total_features_in_transcript}"
            # _temp_solution_right_spacing
        )
        _current += 1

        if run_config.report_transcript_plot_max_samples_per_group is not None:
            limit_transcript_plots_to_samples = flatten_nested_lists([
                sample_names[:run_config.report_transcript_plot_max_samples_per_group]  # Slice each group to max n
                # sample_names
                for sample_names in run_config.sample_groups.values()
            ])
        else:
            limit_transcript_plots_to_samples = None

        plots[f"{fdsa_result.transcript_id}-{fdsa_result.feature_number}"] = {
            "expression": plot_splice_rate(
                fdsa_result,
                norm_gene_counts,
                run_config.sample_groups,
                run_config.group_name_by_sample_name,
                run_config.color_by_group_name,
                run_config.shape_by_group_name,
                show_main_title=False
            ),
            "splice": plot_transcript(
                fdsa_result,
                draw_junctions_with_min_n_occurrences=run_config.report_draw_junctions_min_count,
                show_main_title=False,
                limit_to_samples=limit_transcript_plots_to_samples
            )
        }

    logger.info("... finished")

    logger.info("Generating report...")

    report_html = generate_html_report(
        run_config,
        fdsa_results,
        plots
    )

    output_path = os.path.join(
        # output_dir_absolute, f"Report - {name_output(run_config.report_name, run_config.feature_name)}.html"
        output_dir_absolute, f"{name_output(run_config.report_name, run_config.feature_name)}.html"
    )

    with open(output_path, "w") as output_file:
        output_file.write(report_html)

    logger.info(f"Report written to {output_path}")
    logger.info("...finished")
