
# FEATURECOUNTS DISPATCH AND PARSING GENE COUNT DATA

from dataclasses import dataclass
from typing import Dict, List, Tuple
import os
import subprocess


@dataclass
class FeatureCountsResult:
    """
    Stores results from featureCounts output as:
    1. a dictionary mapping sample names (i.e. column names from header row) to result indices, and
    2. a dictionary mapping gene id to results (counts for each respective sample)
    """
    index_by_sample: Dict[str, int]
    counts_by_gene_id: Dict[str, List[int]]


def get_gene_counts_from_tsv(
    file_path: str,
    gene_id_column_index: int = 0,
    leftmost_counts_column_index: int = 6
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

        header_row = header.replace("\n", "").split("\t")
        index_by_sample_name = {
            sample_name: i for i, sample_name in enumerate(header_row[leftmost_counts_column_index:])
        }

        # Parse remaining lines in file
        counts_by_gene_id = {}
        for line in f:
            row = line.split("\t")
            gene_id = str(row[gene_id_column_index])
            counts_by_gene_id[gene_id] = [int(c) for c in row[leftmost_counts_column_index:]]

    return FeatureCountsResult(index_by_sample_name, counts_by_gene_id)


def run_feature_counts(
    feature_counts_command: str,
    bam_files_dir: str,
    bam_suffix: str,
    reference_gtf_path: str,
    output_path: str,
    paired_end_reads: bool = True,
    threads: int = 1,
    check_exit_code: bool = True
) -> int:

    _suffix_length = len(bam_suffix)
    bam_file_paths = [p for p in os.listdir(bam_files_dir) if p[-_suffix_length:] == bam_suffix]

    paired_flag = "-p " if paired_end_reads else ""

    full_cmd = f"{feature_counts_command} -T {threads} {paired_flag}-t exon -g gene_id --primary " +\
        f"-a {reference_gtf_path} -o {output_path} {' '.join(bam_file_paths)}"

    process = subprocess.run(full_cmd.split(" "), cwd=bam_files_dir, encoding="utf8", check=check_exit_code)

    return process.returncode
