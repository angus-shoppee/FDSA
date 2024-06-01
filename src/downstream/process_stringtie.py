
from typing import Tuple, List, Dict
import csv
import os

from src.utils.general import divide_or_default_zero
# from src.downstream.parse_gtf import STRINGTIE_SOURCE_NAME, SOURCE_INDEX, FEATURE_INDEX, parse_gtf_record
from src.downstream.parse_gtf import get_stringtie_transcripts_from_gtf
import src.config.stringtie_formatted_column_names as cols
import src.config.output_column_names as fase_output_cols


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


def _load_and_sum_results_matrix(
    matrix_path: str
) -> Dict[str, Dict[str, int]]:

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


# def format_stringtie_matrices(
#     gene_counts_path: str,
#     transcript_counts_path: str,
#     merged_gtf_path: str,
#     output_path: str
# ) -> None:
#
#     gene_counts = _load_and_sum_results_matrix(gene_counts_path)
#     transcript_counts = _load_and_sum_results_matrix(transcript_counts_path)
#
#     sample_names = next(iter(gene_counts.values())).keys()
#
#     _output_dir = os.path.dirname(output_path)
#     if not os.path.isdir(_output_dir):
#         os.makedirs(_output_dir)
#
#     with open(output_path, "w") as output_file:
#
#         output_csv = csv.writer(output_file)
#
#         output_csv.writerow(
#             ["seqname", "transcript_id", "gene_id", "ref_gene_id", "ref_gene_name"] +
#             [f"gene.{sample_name}" for sample_name in sample_names] +
#             [f"transcript.{sample_name}" for sample_name in sample_names] +
#             [f"fraction.{sample_name}" for sample_name in sample_names]
#         )
#
#         with open(merged_gtf_path, "r") as gtf_file:
#
#             # _i = 0
#
#             for line in gtf_file:
#
#                 if line[0] == "#":
#                     continue
#
#                 # For a small speed improvement, only parse lines fully if they are StringTie transcripts
#                 _line_split = line.split("\t")
#                 if _line_split[SOURCE_INDEX] == STRINGTIE_SOURCE_NAME and _line_split[FEATURE_INDEX] == "transcript":
#
#                     record = parse_gtf_record(line)
#
#                     gene = gene_counts[record.attributes["gene_id"]]
#                     transcript = transcript_counts[record.attributes["transcript_id"]]
#                     fraction = {
#                         key: divide_or_default_zero(transcript[key], gene[key])
#                         for key in gene.keys()
#                     }
#
#                     output_csv.writerow(
#                         [
#                             record.seqname, record.attributes.get("transcript_id", ""),
#                             record.attributes.get("gene_id", ""), record.attributes.get("ref_gene_id", ""),
#                             record.attributes.get("gene_name", "")
#                             # NOTE: StringTie uses "gene_name" in merged.gtf despite using "ref_gene_name" elsewhere
#                         ] +
#                         [gene[sample_name] for sample_name in sample_names] +
#                         [transcript[sample_name] for sample_name in sample_names] +
#                         [fraction[sample_name] for sample_name in sample_names]
#                     )
#
#                     # print("TRANSCRIPT")
#                     # print(record)
#                     # print("Gene counts:", gene)
#                     # print("Transcript counts:", transcript)
#                     # print("Transcript fractions:", fraction)
#                     # print()
#
#                 # _i += 1
#                 # if _i > 200:
#                 #     break


def format_stringtie_matrices(
    gene_counts_path: str,
    transcript_counts_path: str,
    merged_gtf_path: str,
    output_path: str
) -> None:

    gene_counts = _load_and_sum_results_matrix(gene_counts_path)
    transcript_counts = _load_and_sum_results_matrix(transcript_counts_path)

    sample_names = next(iter(gene_counts.values())).keys()

    _output_dir = os.path.dirname(output_path)
    if not os.path.isdir(_output_dir):
        os.makedirs(_output_dir)

    # print("Parsing transcript details from merged GTF...")

    stringtie_transcripts = get_stringtie_transcripts_from_gtf(merged_gtf_path)

    # print("... done")

    # print("Combining stringtie results...")

    with open(output_path, "w") as output_file:

        output_csv = csv.writer(output_file)

        output_csv.writerow(
            [
                cols.SEQNAME, cols.STRAND, cols.TRANSCRIPT_ID, cols.GENE_ID, cols.REF_GENE_ID, cols.REF_GENE_NAME,
                cols.EXONS
            ] +
            [f"{cols.QUANT_PREFIX_GENE}{sample_name}" for sample_name in sample_names] +
            [f"{cols.QUANT_PREFIX_TRANSCRIPT}{sample_name}" for sample_name in sample_names] +
            [f"{cols.QUANT_PREFIX_FRACTION}{sample_name}" for sample_name in sample_names]
        )

        for gene_id, transcript_record_by_transcript_id in stringtie_transcripts.items():

            for transcript_id, transcript_record in transcript_record_by_transcript_id.items():

                exons_dump = " ".join([
                    f"{start}-{end}"
                    for start, end in transcript_record.exons.values()
                ])

                quant_gene = gene_counts[gene_id]
                quant_transcript = transcript_counts[transcript_id]
                quant_fraction = {
                    key: divide_or_default_zero(quant_transcript[key], quant_gene_value)
                    for key, quant_gene_value in quant_gene.items()
                }

                output_csv.writerow(
                    [
                        transcript_record.seqname, transcript_record.strand, transcript_id, gene_id,
                        transcript_record.attributes.get("ref_gene_id", ""),
                        transcript_record.attributes.get("gene_name", ""),
                        # NOTE: StringTie uses "gene_name" in merged.gtf despite using "ref_gene_name" elsewhere
                        exons_dump
                    ] +
                    [quant_gene[sample_name] for sample_name in sample_names] +
                    [quant_transcript[sample_name] for sample_name in sample_names] +
                    [quant_fraction[sample_name] for sample_name in sample_names]
                )

    # print("... done")


def _get_start_and_end_from_locus(
    locus: str
) -> Tuple[int, int]:

    positions = locus.split("-")

    if len(positions) == 2:
        return int(positions[0]), int(positions[1])

    elif len(positions) == 1:
        return int(positions[0]), int(positions[0])

    else:
        raise ValueError(f"Unable to parse start and end positions from locus: {locus}")


def _calculate_feature_exon_overlap(
    feature_region: Tuple[int, int],
    exons: List[Tuple[int, int]],
    ref_exons: List[Tuple[int, int]]
) -> Tuple[int, float]:

    f_lower_bound = min(feature_region)
    f_upper_bound = max(feature_region)
    f_range = range(f_lower_bound, f_upper_bound + 1)

    ref_e_ranges = [range(min(ref_exon), max(ref_exon) + 1) for ref_exon in ref_exons]
    ref_e_union_set = {n for r in ref_e_ranges for n in r}

    f_nucleotides_within_ref_exons = ref_e_union_set.intersection(f_range)

    e_ranges = [range(min(exon), max(exon) + 1) for exon in exons]
    e_union_set = {n for r in e_ranges for n in r}

    intersection = e_union_set.intersection(f_nucleotides_within_ref_exons)
    intersection_length = len(intersection)

    # overlap_bases, overlap_fraction
    return (
        intersection_length,
        divide_or_default_zero(intersection_length, len(f_nucleotides_within_ref_exons))
    )


def annotate_formatted_stringtie_results(
    formatted_stringtie_output_path: str,
    fase_output_path: str,
    assign_reference_gene: bool = True
) -> None:

    # print("Annotating stringtie results with feature overlap...")

    with open(fase_output_path, "r") as fase_output_file:

        fase_output_csv = csv.reader(fase_output_file)

        header = next(fase_output_csv)
        fase_gene_id_index = header.index(fase_output_cols.GENE_ID)
        fase_ref_exons_index = header.index(fase_output_cols.EXON_POSITIONS)
        fase_feature_region_index = header.index(fase_output_cols.FEATURE_REGION)

        feature_regions: Dict[str, Tuple[int, int]] = {}
        ref_exons: Dict[str, List[Tuple[int, int]]] = {}
        for row in fase_output_csv:
            fase_gene_id = str(row[fase_gene_id_index])
            feature_regions[fase_gene_id] = _get_start_and_end_from_locus(
                row[fase_feature_region_index]
            )
            ref_exons[fase_gene_id] = [
                _get_start_and_end_from_locus(locus)
                for locus in row[fase_ref_exons_index].split(" ")
            ]

    with open(formatted_stringtie_output_path, "r") as stringtie_file:

        stringtie_csv = csv.reader(stringtie_file)

        stringtie_header = next(stringtie_csv)

        stringtie_gene_id_index = stringtie_header.index(cols.GENE_ID)
        stringtie_ref_gene_id_index = stringtie_header.index(cols.REF_GENE_ID)
        stringtie_ref_gene_name_index = stringtie_header.index(cols.REF_GENE_NAME)
        stringtie_exons_index = stringtie_header.index(cols.EXONS)

        stringtie_body = [row for row in stringtie_csv]

    ref_gene_id_by_stringtie_gene_id: Dict[str, str] = {}
    ref_gene_name_by_ref_gene_id: Dict[str, str] = {}
    if assign_reference_gene:

        # Only store 1:1 mappings (store "" in cases where a stringtie gene has >1 corresponding ref gene IDs)
        for row in stringtie_body:

            stringtie_gene_id = row[stringtie_gene_id_index]
            ref_gene_id = row[stringtie_ref_gene_id_index]
            ref_gene_name = row[stringtie_ref_gene_name_index]

            stored_ref_gene_id = ref_gene_id_by_stringtie_gene_id.get(stringtie_gene_id)
            if stored_ref_gene_id is None:
                if ref_gene_id == "":
                    continue  # Do not store blanks if they are the first to be encountered
                ref_gene_id_by_stringtie_gene_id[stringtie_gene_id] = ref_gene_id
            elif stored_ref_gene_id == "":
                pass
            else:
                if ref_gene_id != "" and stored_ref_gene_id != ref_gene_id:
                    ref_gene_id_by_stringtie_gene_id[stringtie_gene_id] = ""

            stored_ref_gene_name = ref_gene_name_by_ref_gene_id.get(ref_gene_id)
            if stored_ref_gene_name is None:
                if ref_gene_name == "":
                    continue   # Do not store blanks if they are the first to be encountered
                if ref_gene_id != "":
                    ref_gene_name_by_ref_gene_id[ref_gene_id] = ref_gene_name
            elif stored_ref_gene_name == "":
                pass
            else:
                if ref_gene_name != "" and ref_gene_id != "" and stored_ref_gene_name != ref_gene_name:
                    ref_gene_name_by_ref_gene_id[ref_gene_id] = ref_gene_name

    with open(formatted_stringtie_output_path, "w") as overwrite_stringtie_file:

        overwrite_stringtie_writer = csv.writer(overwrite_stringtie_file)

        annotated_header = stringtie_header[:stringtie_exons_index + 1] + [
            cols.ANNOTATION_CONTAINS_FEATURE_BASES, cols.ANNOTATION_CONTAINS_FEATURE_FRACTION
        ] + stringtie_header[stringtie_exons_index + 1:]

        overwrite_stringtie_writer.writerow(annotated_header)

        for row in stringtie_body:

            gene_id = row[stringtie_gene_id_index]
            ref_gene_id = row[stringtie_ref_gene_id_index]
            ref_gene_name = row[stringtie_ref_gene_name_index]
            if assign_reference_gene:
                # Use ref gene ID & name from other transcript variants, but only if there is a singular ref gene ID
                if ref_gene_id == "":
                    ref_gene_id = ref_gene_id_by_stringtie_gene_id.get(gene_id, "")
                if ref_gene_name == "":
                    ref_gene_name = ref_gene_name_by_ref_gene_id.get(ref_gene_id, "")
                row[stringtie_ref_gene_id_index] = ref_gene_id
                row[stringtie_ref_gene_name_index] = ref_gene_name

            feature_region = feature_regions.get(ref_gene_id)

            if feature_region is not None:
                overlap_bases, overlap_fraction = _calculate_feature_exon_overlap(
                    feature_region,
                    [_get_start_and_end_from_locus(locus) for locus in row[stringtie_exons_index].split(" ")],
                    ref_exons[ref_gene_id]  # Assume ref_gene_id is a valid key since feature_regions uses same key
                )
                # threshold_pass = 1 if feature_overlap_fraction >= overlap_threshold else 0
            else:
                overlap_bases = ""
                overlap_fraction = ""

            overwrite_stringtie_writer.writerow(
                row[:stringtie_exons_index + 1] +
                [str(overlap_bases), str(overlap_fraction)] +
                row[stringtie_exons_index + 1:]
            )

    # print("... done")
