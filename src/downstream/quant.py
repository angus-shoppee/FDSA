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


from typing import List, Dict, Optional
import logging
import os
import subprocess
import csv
import pandas as pd

from config.quant_merged_column_names import (
    QUANT_MERGED_PREFIX_NORM_CPM, QUANT_MERGED_PREFIX_TRANSCRIPTS_WITHOUT_FEATURE
)
from reporting.process_results import FdsaResult
from downstream.process_stringtie import format_stringtie_matrices, annotate_formatted_stringtie_results


# TODO: Enable automatic deletion of generated files


STRINGTIE_OUTPUT_DIR_NAME = "STRINGTIE"

COMBINED_RESULTS_FILE_NAME = "combined_stringtie_results"


logger = logging.getLogger(__name__)


def quantify_isoforms(
    stringtie_executable_path: str,
    prep_de_script_path: str,
    reference_gtf_path: str,
    filtered_bam_dir: str,
    output_base_dir: str,
    bam_file_ending: str = ".bam",
    read_length: Optional[int] = None,
    fdsa_results_path: Optional[str] = None,
    assign_reference_gene: bool = True,
    threads: int = 1,
    check_exit_code: bool = False
) -> None:

    # TODO: Add option to delete intermediary files & directories upon completion

    if not os.path.isfile(reference_gtf_path):
        raise ValueError(f"Invalid path supplied for reference genome: {reference_gtf_path}")

    if not os.path.isdir(filtered_bam_dir):
        raise ValueError(f"Invalid path supplied for filtered BAM directory: {filtered_bam_dir}")

    stringtie_out_dir = os.path.join(output_base_dir, STRINGTIE_OUTPUT_DIR_NAME)
    if not os.path.isdir(stringtie_out_dir):
        if os.path.exists(stringtie_out_dir):
            raise FileExistsError(f"Could not create output directory ({stringtie_out_dir}): already exists as file")
        os.makedirs(stringtie_out_dir)

    _bam_file_ending_length = len(bam_file_ending)
    bam_abs_paths = [
        os.path.join(os.path.abspath(filtered_bam_dir), f)
        for f in os.listdir(filtered_bam_dir)
        if f[-_bam_file_ending_length:] == bam_file_ending
    ]
    if len(bam_abs_paths) == 0:
        raise ValueError(f"No BAM files found in supplied input directory: {filtered_bam_dir}")

    # Run stringtie individually for each sample

    logger.info("[Stringtie] Assembling transcripts from individual BAM files...")

    assembly_dir = os.path.join(stringtie_out_dir, "assembly")
    stringtie_assembly_output_locations: Dict[str, str] = {}

    # _n = 0
    for bam_path in bam_abs_paths:

        logger.info(os.path.basename(bam_path))

        _basename = os.path.basename(bam_path)
        _output_local_filename = _basename[:_basename.index(".")] + ".gtf"

        stringtie_assembly_output_locations[bam_path] = os.path.join(assembly_dir, _output_local_filename)

        # stringtie -G GENOME -o INDIVIDUAL.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-G", reference_gtf_path,
                "-o", stringtie_assembly_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=stringtie_out_dir,
            encoding="utf8",
            check=check_exit_code
        )

    logger.info("[Stringtie] ... done")

    # Run stringtie in --merge mode to get the combined transcript GTF

    logger.info("[Stringtie] Merging results...")

    # Write merge list for --merge mode
    merge_list_path = os.path.join(assembly_dir, "merge_list.txt")
    with open(merge_list_path, "w") as f:
        f.write("\n".join(list(stringtie_assembly_output_locations.values())))

    merged_gtf_path = os.path.join(stringtie_out_dir, "merged.gtf")

    # stringtie --merge -G GENOME -o MERGED.gtf MERGE_LIST
    subprocess.run(
        [
            stringtie_executable_path,
            "--merge",
            "-G", reference_gtf_path,
            "-o", merged_gtf_path,
            merge_list_path
        ],
        cwd=stringtie_out_dir,
        encoding="utf8",
        check=check_exit_code
    )

    logger.info("[Stringtie] ... done")

    # Run stringtie in -e mode
    # In this step, the combined GTF from stringtie --merge is used with -G rather than the reference genome

    logger.info("[Stringtie] Quantifying transcripts...")

    quantified_dir = os.path.join(stringtie_out_dir, "quantified")
    stringtie_quantified_output_locations: Dict[str, str] = {}

    # _n = 0
    for bam_path in bam_abs_paths:

        _basename = os.path.basename(bam_path)
        _basename_stripped = _basename[:_basename.index(".")]
        _output_local_filename = _basename_stripped + ".quant.gtf"

        stringtie_quantified_output_locations[bam_path] = os.path.join(
            quantified_dir, _basename_stripped, _output_local_filename
        )

        logger.info(_basename)

        # stringtie -e -G MERGED.gtf -o OUT_QUANT.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-e",
                "-G", merged_gtf_path,
                "-o", stringtie_quantified_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=stringtie_out_dir,
            encoding="utf8",
            check=check_exit_code
        )

    logger.info("[Stringtie] ... done")

    # Run prepDE script

    logger.info("[Stringtie] Generating transcript count matrix (prepDE.py3)...")

    prep_de_out_dir = os.path.join(stringtie_out_dir, "prepDE")

    if not os.path.isdir(prep_de_out_dir):
        if os.path.isfile(prep_de_out_dir):
            raise ValueError(f"Could not create prepDE output directory {prep_de_out_dir} - File already exists")
        os.makedirs(prep_de_out_dir)

    prep_de_gene_counts_path = os.path.join(prep_de_out_dir, "gene_count_matrix.csv")
    prep_de_transcript_counts_path = os.path.join(prep_de_out_dir, "transcript_count_matrix.csv")

    # prepDE.py3 -i INPUT [-g GENE_COUNTS.csv] [-t TRANSCRIPT_COUNTS.csv] [-l READ_LENGTH]
    prep_de_cmd_split = [
        "python3", prep_de_script_path,
        "-i", quantified_dir,
        "-g", prep_de_gene_counts_path,
        "-t", prep_de_transcript_counts_path
    ]
    if read_length is not None:
        prep_de_cmd_split += ["-l", str(read_length)]
    subprocess.run(
        prep_de_cmd_split,
        cwd=prep_de_out_dir,
        encoding="utf8",
        check=check_exit_code
    )

    logger.info("[Stringtie] ... done")

    logger.info("[FDSA] Combining stringtie results...")

    formatted_stringtie_output_path = os.path.join(stringtie_out_dir, f"{COMBINED_RESULTS_FILE_NAME}.csv")

    format_stringtie_matrices(
        prep_de_gene_counts_path,
        prep_de_transcript_counts_path,
        merged_gtf_path,
        formatted_stringtie_output_path
    )

    logger.info("[FDSA] ... done")

    if fdsa_results_path is not None:

        logger.info("[FDSA] Annotating stringtie results with feature overlap...")

        annotate_formatted_stringtie_results(
            formatted_stringtie_output_path,
            fdsa_results_path,
            assign_reference_gene=assign_reference_gene
        )

        logger.info("[FDSA] ... done")

    else:

        logger.info("[FDSA] No FDSA output supplied - stringtie results will not be annotated.")


def write_merged_frequencies_and_gene_counts(
    fdsa_results: List[FdsaResult],
    norm_gene_counts: pd.DataFrame,
    output_path: str
) -> None:

    sample_names = list(fdsa_results[0].stringtie_frequencies.keys())

    with open(output_path, "w") as f:

        dump_output_writer = csv.writer(f)

        dump_output_writer.writerow(
            FdsaResult.get_serialized_header_for_info_cols() +
            [f"{QUANT_MERGED_PREFIX_NORM_CPM}.{sample_name}" for sample_name in sample_names] +
            [f"{QUANT_MERGED_PREFIX_TRANSCRIPTS_WITHOUT_FEATURE}.{sample_name}" for sample_name in sample_names]
        )

        for fdsa_result in fdsa_results:

            cpm = norm_gene_counts.loc[fdsa_result.gene_id]

            # stringtie_frequencies_list = [self.stringtie_frequencies[sample_name] for sample_name in sample_names]
            # *[f"{fq:.4f}" for fq in stringtie_frequencies_list]

            dump_output_writer.writerow(
                fdsa_result.serialize_info_cols() +
                # [f"{cpm[sample_name]:.4f}" for sample_name in sample_names] +
                # [f"{fdsa_result.frequencies[sample_name]:.4f}" for sample_name in sample_names]
                [cpm[sample_name] for sample_name in sample_names] +
                [fdsa_result.stringtie_frequencies[sample_name] for sample_name in sample_names]
            )
