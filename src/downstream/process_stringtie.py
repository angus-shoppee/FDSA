# Feature-directed Analysis of Splice Events (FASE)
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


from typing import Tuple, List, Dict
import csv
import os

from src.utils.general import divide_or_default_zero
from src.downstream.parse_gtf import get_stringtie_transcripts_from_gtf, GtfTranscript
import src.config.stringtie_formatted_column_names as cols
import src.config.output_column_names as fase_output_cols


def _load_and_sum_results_matrix(
    matrix_path: str
) -> Dict[str, Dict[str, int]]:

    # Sums counts from multiple lines with the same ID. For example, the following two lines:
    # MSTRG.1,           count_1a, count_2a, ...
    # MSTRG.1|GENE_NAME, count_1b, count_2b, ...
    # Will be stored as {"MSTRG.1": {"sample_1": count_1a + count_1b, "sample_2": count_2a + count_2b, ...}}

    counts: Dict[str, Dict[str, int]] = {}

    with open(matrix_path, "r") as f:

        iterator = csv.reader(f)

        header = next(iterator)
        header_length = len(header)
        sample_names = header[1:]

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
) -> None:

    gene_counts = _load_and_sum_results_matrix(gene_counts_path)
    transcript_counts = _load_and_sum_results_matrix(transcript_counts_path)

    sample_names = next(iter(gene_counts.values())).keys()

    _output_dir = os.path.dirname(output_path)
    if not os.path.isdir(_output_dir):
        os.makedirs(_output_dir)

    # print("Parsing transcript details from merged GTF...")

    stringtie_transcripts: Dict[str, Dict[str, GtfTranscript]] = get_stringtie_transcripts_from_gtf(merged_gtf_path)

    # print("... done")

    # print("Combining stringtie results...")

    with open(output_path, "w") as output_file:

        output_csv = csv.writer(output_file)

        # Write header
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

                # Write data
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
                    [f"{quant_fraction[sample_name]:.4f}" for sample_name in sample_names]
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

    feature_regions: Dict[str, Dict[int, Tuple[int, int]]] = {}  # {Gene ID: {Feature number: (start, end)}}
    ref_exons: Dict[str, List[Tuple[int, int]]] = {}  # {Gene ID: [(start, end), ...]}

    # TODO: Explicitly check that fase_output_path is valid

    with open(fase_output_path, "r") as fase_output_file:

        fase_output_csv = csv.reader(fase_output_file)

        header = next(fase_output_csv)
        fase_gene_id_index = header.index(fase_output_cols.GENE_ID)
        fase_ref_exons_index = header.index(fase_output_cols.EXON_POSITIONS)
        fase_feature_number_index = header.index(fase_output_cols.FEATURE_NUMBER)
        fase_feature_region_index = header.index(fase_output_cols.FEATURE_REGION)

        for row in fase_output_csv:

            fase_gene_id = str(row[fase_gene_id_index])

            ref_exons[fase_gene_id] = [
                _get_start_and_end_from_locus(locus)
                for locus in row[fase_ref_exons_index].split(" ")
            ]

            feature_number = int(row[fase_feature_number_index])
            stored_feature_regions = feature_regions.get(fase_gene_id, dict())
            stored_feature_regions[feature_number] = _get_start_and_end_from_locus(
                row[fase_feature_region_index]
            )
            feature_regions[fase_gene_id] = stored_feature_regions

    with open(formatted_stringtie_output_path, "r") as stringtie_file:

        print(f"Loading {formatted_stringtie_output_path}")

        stringtie_csv = csv.reader(stringtie_file)

        stringtie_header = next(stringtie_csv)

        stringtie_gene_id_index = stringtie_header.index(cols.GENE_ID)
        stringtie_ref_gene_id_index = stringtie_header.index(cols.REF_GENE_ID)
        stringtie_ref_gene_name_index = stringtie_header.index(cols.REF_GENE_NAME)
        stringtie_exons_index = stringtie_header.index(cols.EXONS)

        stringtie_body = [row for row in stringtie_csv]

        print("Finished loading")

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

    # Write to file
    with open(formatted_stringtie_output_path, "w") as overwrite_stringtie_file:

        overwrite_stringtie_writer = csv.writer(overwrite_stringtie_file)

        # Write header
        overwrite_stringtie_writer.writerow(
            stringtie_header[:stringtie_exons_index + 1] +
            [
                cols.ANNOTATION_FEATURE_REGIONS,
                cols.ANNOTATION_CONTAINS_FEATURE_BASES,
                cols.ANNOTATION_CONTAINS_FEATURE_FRACTION
            ] +
            stringtie_header[stringtie_exons_index + 1:]
        )

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

            feature_region_by_feature_number = feature_regions.get(ref_gene_id, None)
            feature_regions_string_dump = ""
            overlap_bases_string_dump = ""
            overlap_fraction_string_dump = ""

            if feature_region_by_feature_number is not None:

                for feature_number, feature_region in feature_region_by_feature_number.items():

                    feature_regions_string_dump += f"f{feature_number}:{feature_region[0]}-{feature_region[1]} "

                    overlap_bases, overlap_fraction = _calculate_feature_exon_overlap(
                        feature_region,
                        [_get_start_and_end_from_locus(locus) for locus in row[stringtie_exons_index].split(" ")],
                        ref_exons[ref_gene_id]  # Assume ref_gene_id is a valid key since feature_regions uses same key
                    )
                    overlap_bases_string_dump += f"f{feature_number}:{overlap_bases} "
                    overlap_fraction_string_dump += f"f{feature_number}:{overlap_fraction:.2f} "

            # Write data
            overwrite_stringtie_writer.writerow(
                row[:stringtie_exons_index + 1] +
                [
                    feature_regions_string_dump.rstrip(),
                    overlap_bases_string_dump.rstrip(),
                    overlap_fraction_string_dump.rstrip()
                ] +
                row[stringtie_exons_index + 1:]
            )

    # print("... done")
