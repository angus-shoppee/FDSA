
from typing import Dict
import csv
import os

from src.utils.general import divide_or_default_zero
from src.downstream.parse_gtf import STRINGTIE_SOURCE_NAME, SOURCE_INDEX, FEATURE_INDEX, parse_gtf_record


# def _load_results_matrix(matrix_path: str) -> Dict[str: Dict[str, int]]:
#
#     with open(matrix_path, "r") as f:
#
#         iterator = csv.reader(matrix_path)
#
#         header = next(iterator)
#         header_length = len(header)
#
#         return {
#             row[0]: {
#                 header[i]: row[i]
#                 for i in range(1, header_length)
#             }
#             for row in iterator
#         }


def _load_and_sum_results_matrix(matrix_path: str) -> Dict[str, Dict[str, int]]:

    # Sums counts from multiple lines with the same ID. For example, the following two lines:
    # MSTRG.1,           count_1a, count_2a, ...
    # MSTRG.1|GENE_NAME, count_1b, count_2b, ...
    # Would be stored as {"MSTRG.1": {"sample_1": count_1a + count_1b, "sample_2": count_2a + count_2b, ...}}

    with open(matrix_path, "r") as f:

        iterator = csv.reader(f)

        header = next(iterator)
        header_length = len(header)
        sample_names = header[1:]

        counts: Dict[str, Dict[str, int]] = {}

        for row in iterator:

            identifier = row[0].partition("|")[0]

            if counts.get(identifier) is None:

                counts[identifier] = {
                    header[i]: int(row[i])
                    for i in range(1, header_length)
                }

            else:

                stored = counts[identifier]

                for i, sample_name in enumerate(sample_names):

                    stored[sample_name] += int(row[i + 1])  # Add 1 to index to account for identifier column

        return counts


def format_stringtie_matrices(
    gene_counts_path: str,
    transcript_counts_path: str,
    merged_gtf_path: str,
    output_path: str
):

    gene_counts = _load_and_sum_results_matrix(gene_counts_path)
    transcript_counts = _load_and_sum_results_matrix(transcript_counts_path)

    sample_names = next(iter(gene_counts.values())).keys()

    _output_dir = os.path.dirname(output_path)
    if not os.path.isdir(_output_dir):
        os.makedirs(_output_dir)

    with open(output_path, "w") as output_file:

        output_csv = csv.writer(output_file)

        output_csv.writerow(
            ["seqname", "transcript_id", "gene_id", "ref_gene_id", "ref_gene_name"] +
            [f"gene.{sample_name}" for sample_name in sample_names] +
            [f"transcript.{sample_name}" for sample_name in sample_names] +
            [f"fraction.{sample_name}" for sample_name in sample_names]
        )

        with open(merged_gtf_path, "r") as gtf_file:

            # _i = 0

            for line in gtf_file:

                if line[0] == "#":
                    continue

                # For a small speed improvement, only parse lines fully if they are the desired feature type
                _line_split = line.split("\t")
                if _line_split[SOURCE_INDEX] == STRINGTIE_SOURCE_NAME and _line_split[FEATURE_INDEX] == "transcript":

                    record = parse_gtf_record(line)

                    gene = gene_counts[record.attributes["gene_id"]]
                    transcript = transcript_counts[record.attributes["transcript_id"]]
                    fraction = {
                        key: divide_or_default_zero(transcript[key], gene[key])
                        for key in gene.keys()
                    }

                    output_csv.writerow(
                        [
                            record.seqname, record.attributes.get("transcript_id", ""),
                            record.attributes.get("gene_id", ""), record.attributes.get("ref_gene_id", ""),
                            record.attributes.get("gene_name", "")
                            # NOTE: StringTie uses "gene_name" in merged.gtf despite using "ref_gene_name" elsewhere
                        ] +
                        [gene[sample_name] for sample_name in sample_names] +
                        [transcript[sample_name] for sample_name in sample_names] +
                        [fraction[sample_name] for sample_name in sample_names]
                    )

                    # if record.attributes.get("gene_id") == "MSTRG.620":
                    print("TRANSCRIPT")
                    print(record)
                    print("Gene counts:", gene)
                    print("Transcript counts:", transcript)
                    print("Transcript fractions:", fraction)
                    print()

                # _i += 1
                # if _i > 200:
                #     break
