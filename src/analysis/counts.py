
# FEATURECOUNTS DISPATCH AND PARSING GENE COUNT DATA

from dataclasses import dataclass
from typing import Dict, List
import os
import subprocess
import pandas as pd
import conorm


DEFAULT_FEATURE_COUNTS_OUTPUT_GENE_ID_COLUMN_INDEX = 0
DEFAULT_FEATURE_COUNTS_OUTPUT_LEFTMOST_COUNTS_COLUMN_INDEX = 6


@dataclass
class FeatureCountsResult:
    """
    Stores results from featureCounts output as:
    1. a dictionary mapping sample names (i.e. column names from header row) to result indices, and
    2. a dictionary mapping gene id to results (counts for each respective sample)
    """
    index_by_sample: Dict[str, int]
    counts_by_gene_id: Dict[str, List[int]]


def get_tmm_cpm_from_gene_counts(
    feature_counts_result: FeatureCountsResult
) -> pd.DataFrame:

    sample_names_ordered = sorted(
        [(i, sample_name) for (sample_name, i) in feature_counts_result.index_by_sample.items()],
        key=lambda x: x[0]
    )

    row_names, rows = [], []
    for gene_id, counts in feature_counts_result.counts_by_gene_id.items():
        row_names.append(gene_id)
        rows.append(counts)

    counts = pd.DataFrame(
        rows,
        index=row_names,
        columns=[x[1] for x in sample_names_ordered]
    )

    norm_factors = conorm.tmm_norm_factors(counts)

    tmm_cpm = conorm.cpm(counts, norm_factors=norm_factors)

    return tmm_cpm


def get_gene_counts_from_tsv(
    file_path: str,
    bam_ending: str = "",
    gene_id_column_index: int = DEFAULT_FEATURE_COUNTS_OUTPUT_GENE_ID_COLUMN_INDEX,
    leftmost_counts_column_index: int = DEFAULT_FEATURE_COUNTS_OUTPUT_LEFTMOST_COUNTS_COLUMN_INDEX
) -> FeatureCountsResult:

    """Parses a featureCounts output file and returns data as a FeatureCountsResult instance"""

    # Load counts file
    with open(file_path, "r") as f:

        # Check for metadata row
        line_0 = f.readline()
        if line_0[0] == "#":
            header = f.readline()  # Progress to next row if metadata row is present
        else:
            header = line_0  # Otherwise, the first line of the file is the header row

        if bam_ending:
            header = header.replace(bam_ending, "")

        header_row = header.replace("\n", "").split("\t")

        index_by_sample_name = {
            sample_name: i for i, sample_name in enumerate(header_row[leftmost_counts_column_index:])
        }

        # Parse remaining lines in file
        counts_by_gene_id = {}
        for line in f:
            row = line.replace("\n", "").split("\t")
            gene_id = str(row[gene_id_column_index])
            counts_by_gene_id[gene_id] = [int(c) for c in row[leftmost_counts_column_index:]]

    return FeatureCountsResult(index_by_sample_name, counts_by_gene_id)


def run_feature_counts(
    feature_counts_executable: str,
    bam_files_dir: str,
    bam_ending: str,
    reference_gtf_path: str,
    output_path: str,
    paired_end_reads: bool = True,
    primary_alignment_only: bool = False,
    threads: int = 1,
    check_exit_code: bool = True
) -> int:

    _ending_length = len(bam_ending)
    bam_file_paths = [p for p in os.listdir(bam_files_dir) if p[-_ending_length:] == bam_ending]

    paired_flag = "-p " if paired_end_reads else ""
    primary_flag = "--primary " if primary_alignment_only else ""

    full_cmd = f"{feature_counts_executable} -T {threads} {paired_flag}-t exon -g gene_id {primary_flag}" +\
        f"-a {reference_gtf_path} -o {output_path} {' '.join(bam_file_paths)}"

    process = subprocess.run(full_cmd.split(" "), cwd=bam_files_dir, encoding="utf8", check=check_exit_code)

    return process.returncode
