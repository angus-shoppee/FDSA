
# FASE CONFIG & CORE FUNCTIONS

from dataclasses import dataclass
from typing import Tuple, List, Dict, Any
from multiprocessing import Pool
import os
import csv
import configparser
import time  # DEBUG

from transcript import TranscriptRecord, TranscriptLibrary
from experiment import Sample
from counts import FeatureCountsResult
from splice import get_splice_junctions_from_sample


class FaseConfig:
    allowed_species: List[str]
    biomart_name_for_species: Dict[str, str]

    def __init__(self, internal_config_path):

        config = configparser.ConfigParser()
        config.read(internal_config_path)

        self.allowed_species = config["DEFAULT"]["AllowedSpecies"].split(" ")
        self.biomart_name_for_species = {k: v for k, v in [
            s.split(":") for s in config["DEFAULT"]["BiomartNameForSpecies"].split(" ")
        ]}


def sanitize_string_for_filename(s: str) -> str:
    return "".join(c for c in s if c.isalnum() or c in (".", "_", "-", " ")).rstrip()


def calculate_optc(n_occurrences: int, n_reads: int) -> float:

    return 0.0 if not all((n_occurrences, n_reads)) else float((n_occurrences * 1000) / n_reads)


@dataclass
class ExonicOverlap:
    bases: int
    fraction: float


def calculate_exonic_overlap(
    transcript: TranscriptRecord,
    f_pos: Tuple[int, int],
    j_pos: Tuple[int, int]
) -> ExonicOverlap:

    f_start, f_end = f_pos
    j_start, j_end = j_pos

    j_range = range(j_start, j_end + 1)
    f_range = range(f_start, f_end + 1) if f_start <= f_end else range(f_end, f_start + 1)  # Account for reverse strand

    e_ranges = [range(e.start, e.end + 1) for e in transcript.exons]

    # print(e_ranges)

    # Ugly way to make a set from a list of ranges... is there a more pythonic way to do this?
    e_union_set = set()
    for r in e_ranges:
        e_union_set = e_union_set.union(set(r))

    e_j_intersect = e_union_set.intersection(j_range)

    e_f_intersect = e_union_set.intersection(f_range)
    length_of_feature = len(e_f_intersect)

    exonic_overlap = e_j_intersect.intersection(e_f_intersect)
    length_of_overlap = len(exonic_overlap)

    overlap_fraction = 0 if not length_of_overlap else length_of_overlap / length_of_feature

    return ExonicOverlap(length_of_overlap, overlap_fraction)


@dataclass
class OverlappingJunctionsInfo:
    sample_name: str
    number_occurrences: int
    frequency_occurrences: float
    junctions_loci: str


def get_splice_analysis_result_for_sample(
    sample: Sample,
    transcript: TranscriptRecord,
    gene_counts: FeatureCountsResult,
    overlap_threshold: float,
    f_start_genomic: int,
    f_end_genomic: int,
    primary_alignment_only: bool,
    experimental_allow_junction_start_outside_transcript: bool = False,
    experimental_allow_junction_end_outside_transcript: bool = False
) -> OverlappingJunctionsInfo:

    number_occurrences = 0

    overlapping_j_buffer = f"| {sample.name} "

    # Load junctions in this locus
    # _start_j_ms = time.time_ns() / 1_000_000
    junctions = get_splice_junctions_from_sample(
        sample,
        transcript.seqname,
        transcript.start,
        transcript.end,
        primary_alignment_only=primary_alignment_only
    )
    # _finish_j_ms = time.time_ns() / 1_000_000
    # print(f"INFO: ran get_splice_junctions_from_sam in {_finish_j_ms - _start_j_ms} ms")

    for junction in junctions:

        # Test junction bounds
        if junction.start < transcript.start or junction.start > transcript.end:
            if not experimental_allow_junction_start_outside_transcript:
                continue
        if junction.end < transcript.start or junction.end > transcript.end:
            if not experimental_allow_junction_end_outside_transcript:
                continue
        if junction.start <= transcript.start and junction.end >= transcript.end:
            continue

        # Determine overlap between the junction and the locus of the feature
        overlap = calculate_exonic_overlap(
            transcript,
            (f_start_genomic, f_end_genomic),
            (junction.start, junction.end)
        )

        # Add n unique to tally of occurrences if overlap exceeds user-defined threshold
        if overlap.fraction >= overlap_threshold:
            overlapping_j_buffer += f"{junction.n_unique}*[{junction.start}-{junction.end}] "

            number_occurrences += junction.n_unique

    # Set OPTC after tallying occurrences in all relevant junctions
    frequency_occurrences = calculate_optc(
        number_occurrences,
        gene_counts.counts_by_gene_id[
            transcript.gene_id
        ][gene_counts.index_by_sample[f"{sample.name}{sample.suffix}"]]
    )

    return OverlappingJunctionsInfo(
        sample.name,
        number_occurrences,
        frequency_occurrences,
        overlapping_j_buffer
    )


# def _print_transcript_debug_info(
#     transcript: TranscriptRecord,
#     feature_substring: str,
#     splice_junction_loci: str
# ) -> None:
#
#     print("-" * 40)
#
#     print(f"{transcript.transcript_id} / {transcript.gene_name}")
#     print(f"Strand: ({transcript.strand})")
#     print(f"Transcript exons: {[(e.start, e.end) for e in transcript.exons]}")
#     print(f"GBSeq exons: {[(e.start, e.end) for e in transcript.gbseq.exons]}")
#     print(f"Exon map: {transcript.exon_map}")
#     print(f"Overlapping junctions: {splice_junction_loci}")
#
#     for feature in transcript.gbseq.features:
#         if feature.has_qual_value_containing(feature_substring):
#             try:
#                 g_start = transcript.genomic_position_from_local(feature.start)
#                 g_end = transcript.genomic_position_from_local(feature.end)
#                 print(f"Feature with '{feature_substring}' qual: Local = " +
#                       f"{str(feature.start) + ', ' + str(feature.end)}")
#                 print(f"Feature with '{feature_substring}' qual: Genomic = " +
#                       f"{str(g_start) + ', ' + str(g_end)}")
#             except ValueError:
#                 print("Unable to convert local positions to global")
#
#     print("-" * 40)


def perform_splice_analysis(
    features_to_analyze: List[str],
    overlap_threshold: float,
    samples: Dict[str, Sample],
    transcript_library: TranscriptLibrary,
    gene_counts: FeatureCountsResult,
    output_dir: str,
    max_n_processes: int = 1,
    ns_threshold_to_use_multiprocessing: int = 100_000_000,
    primary_alignment_only: bool = False,
    verbose: bool = True,
) -> None:

    def print_if_verbose(s: Any):
        if verbose:
            print(s)

    if not 0.0 <= overlap_threshold <= 1.0:
        raise ValueError(f"overlap_threshold must be between 0 and 1. Got: {overlap_threshold}")
    if not os.path.isdir(output_dir):
        raise ValueError(f"output_directory must be a valid directory. Got: {output_dir}")

    _analysis_start_ms = time.time_ns() / 1_000_000

    n_samples = len(samples)

    n_processes = min(len(samples), max_n_processes)

    pool = Pool(processes=n_processes)
    print_if_verbose(f"Created worker pool with {n_processes} processes for parallel BAM file parsing\n")

    for feature_substring in features_to_analyze:

        print_if_verbose(f"Performing splice analysis for term \"{feature_substring}\"...")

        skipped_exon_match_failure = 0

        with open(os.path.join(output_dir, f"{sanitize_string_for_filename(feature_substring)}.csv"), "w") as f:

            sample_names_alphabetical = sorted(samples.keys())
            header = [
                "Transcript ID",
                "Gene name",
                "Transcript genomic start position",
                "Exon positions",
                "N instances of feature within transcript",
                "Feature number within transcript",
                "Feature region",
                "Overlapping junctions in samples"
            ] + sample_names_alphabetical + [f"OPTC {sample_name}" for sample_name in sample_names_alphabetical]

            out_csv = csv.writer(f)
            out_csv.writerow(header)

            for _i, transcript in enumerate(transcript_library.get_all_transcripts()):

                # TODO: Exclude "redundant" transcripts with duplicate or matching annotation.
                #       Skipping prior to performing splice analysis will result in large time savings and additionally
                #       reduce size and redundancy of output.

                # START DEBUG BLOCK
                if transcript.gene_name not in ("Cacna1d", "Cd74", "Gpr6", "Pkd1"):
                    continue
                if _i > 5000:
                    break
                # print("REGION:", f"{transcript.seqname}:{transcript.start}-{transcript.end}")
                # END OF DEBUG BLOCK

                # # START DEBUG BLOCK 2
                # if transcript.gene_name not in ("Cd274", "Tnfrsf1b", "Bcl2"):
                #     continue
                # # END OF DEBUG BLOCK 2

                # # START OF DEBUG BLOCK 3
                # if _i > 1000:
                #     break
                #     # if transcript.gene_name not in ("Cacna1d", "Cd74", "Gpr6", "Pkd1"):
                #     #     continue
                #     # if _i > 5000:
                #     #     break
                # # END OF DEBUG BLOCK 3

                _start_ms = time.time_ns() / 1_000_000

                if verbose:
                    if _i % 1000 == 0:
                        print(f"Progress: {_i}/{transcript_library.number_of_transcripts}")

                analysis_features = transcript.gbseq.get_analysis_features()

                if feature_substring in analysis_features.keys():

                    features = analysis_features[feature_substring]
                    n_features_in_transcript = len(features)

                    for feature_index, feature in enumerate(features):

                        try:
                            f_start_genomic = transcript.genomic_position_from_local(feature.start)
                            f_end_genomic = transcript.genomic_position_from_local(feature.end)
                        except ValueError:
                            skipped_exon_match_failure += 1
                            continue

                        feature_region = f"{transcript.seqname}({transcript.strand}):{f_start_genomic}-{f_end_genomic}"
                        exon_positions = " ".join([f"{e.start}-{e.end}" for e in transcript.exons])
                        overlapping_junctions_loci = ""

                        raw_number_occurrences = {}  # Dict[str, int]
                        frequency_occurrences = {}  # Dict[str, float]

                        samples_alphabetical = [samples[sample_name] for sample_name in sample_names_alphabetical]

                        # Test time taken to parse first BAM file and decide whether parallelization is necessary
                        first_bam_start_ns = time.time_ns()
                        first_result = get_splice_analysis_result_for_sample(
                            samples_alphabetical[0],
                            transcript,
                            gene_counts,
                            overlap_threshold,
                            f_start_genomic,
                            f_end_genomic,
                            primary_alignment_only
                        )
                        first_bam_finish_ns = time.time_ns()

                        # Proceed with either serial function calls or pool.starmap if there is more than one BAM file
                        if n_samples > 1:

                            if first_bam_finish_ns - first_bam_start_ns > ns_threshold_to_use_multiprocessing:
                                args_to_pool = [[
                                    sample,
                                    transcript,
                                    gene_counts,
                                    overlap_threshold,
                                    f_start_genomic,
                                    f_end_genomic,
                                    primary_alignment_only
                                ] for sample in samples_alphabetical]
                                remaining_results = pool.starmap(
                                    get_splice_analysis_result_for_sample,
                                    args_to_pool
                                )

                            else:
                                remaining_results = [
                                    get_splice_analysis_result_for_sample(
                                        sample,
                                        transcript,
                                        gene_counts,
                                        overlap_threshold,
                                        f_start_genomic,
                                        f_end_genomic,
                                        primary_alignment_only
                                    ) for sample in samples_alphabetical[1:]
                                ]

                        else:

                            remaining_results = []

                        results = [first_result] + remaining_results

                        for result in results:
                            raw_number_occurrences[result.sample_name] = result.number_occurrences
                            frequency_occurrences[result.sample_name] = result.frequency_occurrences
                            overlapping_junctions_loci += result.junctions_loci

                        # _print_transcript_debug_info(transcript, feature_substring, overlapping_junctions_loci)

                        # Write to output (one row per feature per transcript)
                        out_csv.writerow(
                            [
                                transcript.transcript_id,
                                transcript.gene_name,
                                transcript.start,
                                exon_positions,
                                n_features_in_transcript,
                                feature_index + 1,
                                feature_region,
                                overlapping_junctions_loci
                            ] + [
                                raw_number_occurrences[sample_name] for sample_name in sample_names_alphabetical
                            ] + [
                                frequency_occurrences[sample_name] for sample_name in sample_names_alphabetical
                            ]
                        )

                _finish_ms = time.time_ns() / 1_000_000
                print(
                    f"INFO: Processed transcript ({transcript.transcript_id} / {transcript.gene_name}) " +
                    f"in {_finish_ms - _start_ms} ms"
                )

        _analysis_finish_ms = time.time_ns() / 1_000_000
        _analysis_runtime_minutes = (_analysis_finish_ms - _analysis_start_ms) / (1000 * 60)

        print(f"TIME: Completed in {_analysis_runtime_minutes} minutes")

        print_if_verbose(
            f"...finished " +
            f"({skipped_exon_match_failure} skipped due to failed matching between reference and annotation exons)\n"
        )
