
# FASE CONFIG & CORE FUNCTIONS

from dataclasses import dataclass
from typing import Tuple, List, Dict, Any, Union
from multiprocessing import Pool
import os
import csv
import configparser
import time  # DEBUG

from transcript import TranscriptRecord, TranscriptLibrary
from experiment import Sample
from counts import FeatureCountsResult
from splice import SpliceJunction, get_splice_junctions_from_sample


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


def set_analysis_features(
    features_to_analyze: List[str],
    annotated_transcript_library: TranscriptLibrary,
    only_use_longest_annotated_transcript: bool = False,
    skip_transcripts_with_redundant_feature_annotation: bool = False
) -> None:

    for feature_substring in features_to_analyze:

        if skip_transcripts_with_redundant_feature_annotation or only_use_longest_annotated_transcript:

            transcripts_by_gene = annotated_transcript_library.get_transcripts_for_all_genes()

            for gene_name in transcripts_by_gene.keys():

                transcripts_ordered_by_length = sorted(
                    transcripts_by_gene[gene_name].values(),
                    key=lambda x: x.get_length(),
                    reverse=True
                )

                if only_use_longest_annotated_transcript:

                    reached_transcript_with_relevant_annotation = False
                    for transcript in transcripts_ordered_by_length:

                        for gbfeature in transcript.gbseq.features:
                            if gbfeature.has_qual_value_containing(feature_substring):
                                reached_transcript_with_relevant_annotation = True
                                transcript.gbseq.set_analysis_feature(feature_substring, gbfeature)

                        if reached_transcript_with_relevant_annotation:
                            break

                else:

                    # Get relevant features and store against transcript_id
                    analysis_features = {
                        transcript.transcript_id: list() for transcript in transcripts_ordered_by_length
                    }
                    for transcript_id, transcript in transcripts_by_gene[gene_name].items():
                        for gbfeature in transcript.gbseq.features:
                            if gbfeature.has_qual_value_containing(feature_substring):
                                analysis_features[transcript_id].append(gbfeature)

                    # Only use set_analysis_feature for transcripts with unique annotation, defaulting to longest
                    encountered = []
                    for transcript in transcripts_ordered_by_length:
                        summary_of_analysis_features = [
                            (
                                gbfeature.start, gbfeature.end
                            ) for gbfeature in analysis_features[transcript.transcript_id]
                        ]
                        if summary_of_analysis_features not in encountered:
                            encountered.append(summary_of_analysis_features)
                            for gbfeature in analysis_features[transcript.transcript_id]:
                                transcript.gbseq.set_analysis_feature(feature_substring, gbfeature)

        else:

            for transcript in annotated_transcript_library.get_all_transcripts():
                for gbfeature in transcript.gbseq.features:
                    if gbfeature.has_qual_value_containing(feature_substring):
                        transcript.gbseq.set_analysis_feature(feature_substring, gbfeature)


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
    junctions: List[SpliceJunction],
    transcript: TranscriptRecord,
    n_gene_counts: int,
    overlap_threshold: float,
    f_start_genomic: int,
    f_end_genomic: int,
    experimental_allow_junction_start_outside_transcript: bool = False,
    experimental_allow_junction_end_outside_transcript: bool = False
) -> OverlappingJunctionsInfo:

    number_occurrences = 0

    overlapping_j_buffer = f"| {sample.name} "

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
        n_gene_counts
        # gene_counts.counts_by_gene_id[
        #     transcript.gene_id
        # ][gene_counts.index_by_sample[f"{sample.name}{sample.suffix}"]]
    )

    return OverlappingJunctionsInfo(
        sample.name,
        number_occurrences,
        frequency_occurrences,
        overlapping_j_buffer
    )


def perform_splice_analysis(
    features_to_analyze: List[str],
    overlap_threshold: float,
    samples: Dict[str, Sample],
    transcript_library: TranscriptLibrary,
    gene_counts: FeatureCountsResult,
    output_dir: str,
    max_n_features_in_transcript: Union[None, int] = None,
    max_n_processes: int = 1,
    min_bam_parse_time_ns_to_use_multiprocessing: int = 100_000_000,
    min_total_n_junctions_to_use_multiprocessing: int = 300,
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

    if n_processes > 1:
        pool = Pool(processes=n_processes)
        print_if_verbose(f"Created worker pool with {n_processes} processes\n")
    else:
        pool = None

    for feature_substring in features_to_analyze:

        # TODO: Refactor to iterate over transcripts & samples, and perform analysis for all feature substrings with
        #       only one BAM file query. i.e.:
        #       get_splice_junctions_from_sample --> call get_splice_analysis_result_for_sample for each feature
        #       (will require refactoring get_splice_analysis_result_for_sample to accept junctions rather than sample)

        print_if_verbose(f"Performing splice analysis for term \"{feature_substring}\"...")

        skipped_exon_match_failure = 0

        with open(os.path.join(output_dir, f"{sanitize_string_for_filename(feature_substring)}.csv"), "w") as f:

            sample_names_alphabetical = sorted(samples.keys())
            samples_alphabetical = [samples[sample_name] for sample_name in sample_names_alphabetical]

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

            _info_time_bam_read = False
            _info_time_transcripts = False

            _progress = 0
            for transcripts_by_id in transcript_library.get_transcripts_for_all_genes().values():

                if not transcripts_by_id:
                    continue

                transcripts = list(transcripts_by_id.values())

                # # START DEBUG BLOCK
                # if transcripts[0].gene_name not in ("Cacna1d", "Cd74", "Gpr6", "Pkd1"):
                #     _progress += 1
                #     continue
                # if _progress > 5000:
                #     break
                # # print("REGION:", f"{transcript.seqname}:{transcript.start}-{transcript.end}")
                # # END OF DEBUG BLOCK

                n_gene_counts_by_sample_name = {
                    sample.name: gene_counts.counts_by_gene_id[
                        transcripts[0].gene_id
                    ][gene_counts.index_by_sample[f"{sample.name}{sample.suffix}"]] for sample in samples_alphabetical
                }

                chromosome = transcripts[0].seqname
                gene_region_start = min([transcript.start for transcript in transcripts])
                gene_region_end = min([transcript.end for transcript in transcripts])

                if _info_time_bam_read:
                    _bam_read_start_ms = time.time_ns() / 1_000_000

                # Test time taken to parse junctions from first BAM file and decide whether parallelization is necessary
                first_bam_start_ns = time.time_ns()
                first_junctions = get_splice_junctions_from_sample(
                    samples_alphabetical[0],
                    chromosome,
                    gene_region_start,
                    gene_region_end,
                    primary_alignment_only
                )
                first_bam_finish_ns = time.time_ns()

                # Proceed with either serial function calls or pool.starmap if there is more than one BAM file
                if n_samples > 1:

                    if all([
                        pool is not None,
                        first_bam_finish_ns - first_bam_start_ns > min_bam_parse_time_ns_to_use_multiprocessing
                    ]):
                        args_to_pool = [[
                            sample,
                            chromosome,
                            gene_region_start,
                            gene_region_end,
                            primary_alignment_only
                        ] for sample in samples_alphabetical]
                        remaining_junctions = pool.starmap(
                            get_splice_junctions_from_sample,
                            args_to_pool
                        )

                    else:
                        remaining_junctions = [
                            get_splice_junctions_from_sample(
                                sample,
                                chromosome,
                                gene_region_start,
                                gene_region_end,
                                primary_alignment_only
                            ) for sample in samples_alphabetical[1:]
                        ]

                else:

                    remaining_junctions = []

                all_junctions = [first_junctions] + remaining_junctions

                junctions_by_sample_name = {
                    sample_names_alphabetical[i]: all_junctions[i] for i in range(len(sample_names_alphabetical))
                }

                if _info_time_bam_read:
                    _bam_read_finish_ms = time.time_ns() / 1_000_000
                    print(
                        f"INFO: parsed splice junctions from BAM files for gene {transcripts[0].gene_name} in " +
                        f"{_bam_read_finish_ms - _bam_read_start_ms} ms"
                    )

                for transcript in transcripts:

                    if verbose:
                        if _progress % 1000 == 0:
                            print(f"Progress: {_progress}/{transcript_library.number_of_transcripts}")

                    analysis_features = transcript.gbseq.get_analysis_features()

                    # REFACTOR --> if any([s in analysis_features.keys() for s in features_to_analyze]):
                    if feature_substring in analysis_features.keys():

                        if _info_time_transcripts:
                            _transcript_start_ms = time.time_ns() / 1_000_000

                        features = analysis_features[feature_substring]
                        n_features_in_transcript = len(features)

                        if max_n_features_in_transcript:
                            if n_features_in_transcript > max_n_features_in_transcript:
                                _progress += 1
                                continue

                        for feature_index, feature in enumerate(features):

                            # TODO: Implement multiprocessing of multiple features at a time for larger worker pools

                            try:
                                f_start_genomic = transcript.genomic_position_from_local(feature.start)
                                f_end_genomic = transcript.genomic_position_from_local(feature.end)
                            except ValueError:
                                skipped_exon_match_failure += 1
                                continue

                            feature_region = \
                                f"{transcript.seqname}({transcript.strand}):{f_start_genomic}-{f_end_genomic}"
                            exon_positions = " ".join([f"{e.start}-{e.end}" for e in transcript.exons])
                            overlapping_junctions_loci = ""

                            raw_number_occurrences = {}  # Dict[str, int]
                            frequency_occurrences = {}  # Dict[str, float]

                            if all([
                                pool is not None,
                                sum(
                                    [len(j) for j in junctions_by_sample_name.values()]
                                ) > min_total_n_junctions_to_use_multiprocessing
                            ]):
                                args_to_pool = [[
                                    sample,
                                    junctions_by_sample_name[sample.name],
                                    transcript,
                                    n_gene_counts_by_sample_name[sample.name],
                                    overlap_threshold,
                                    f_start_genomic,
                                    f_end_genomic
                                ] for sample in samples_alphabetical]
                                results = pool.starmap(
                                    get_splice_analysis_result_for_sample,
                                    args_to_pool
                                )
                            else:
                                results = [get_splice_analysis_result_for_sample(
                                    sample,
                                    junctions_by_sample_name[sample.name],
                                    transcript,
                                    n_gene_counts_by_sample_name[sample.name],
                                    overlap_threshold,
                                    f_start_genomic,
                                    f_end_genomic
                                ) for sample in samples_alphabetical]

                            for result in results:
                                raw_number_occurrences[result.sample_name] = result.number_occurrences
                                frequency_occurrences[result.sample_name] = result.frequency_occurrences
                                overlapping_junctions_loci += result.junctions_loci

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

                        if _info_time_transcripts:
                            _transcript_finish_ms = time.time_ns() / 1_000_000
                            print(
                                f"INFO: Processed transcript ({transcript.transcript_id} / {transcript.gene_name}) " +
                                f"in {_transcript_finish_ms - _transcript_start_ms} ms"
                            )

                    _progress += 1

        _analysis_finish_ms = time.time_ns() / 1_000_000
        _analysis_runtime_minutes = (_analysis_finish_ms - _analysis_start_ms) / (1000 * 60)

        print(f"TIME: Completed in {_analysis_runtime_minutes} minutes")

        print_if_verbose(
            f"...finished " +
            f"({skipped_exon_match_failure} skipped due to failed matching between reference and annotation exons)\n"
        )
