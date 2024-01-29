
# FASE CONFIG & CORE FUNCTIONS

from dataclasses import dataclass
from typing import Tuple, List, Dict, Any, Union
from multiprocessing import Pool
import os
import csv
import time
from statistics import mean

from src.config import output_column_names as cols
from src.analysis.transcript import TranscriptRecord, TranscriptLibrary
from src.analysis.experiment import Sample
from src.analysis.splice import DEFAULT_MAPQ_FOR_UNIQUE_MAPPING, SpliceJunction, get_splice_junctions_from_sample


def sanitize_string_for_filename(s: str) -> str:
    return "".join(c for c in s if c.isalnum() or c in (".", "_", "-", " ")).rstrip()


def name_output_file(
    run_name: str,
    feature_name: str
) -> str:
    return sanitize_string_for_filename(f"{run_name} - {feature_name}")


def set_analysis_features(
    features_to_analyze: List[str],
    annotated_transcript_library: TranscriptLibrary,
    only_use_longest_annotated_transcript: bool = True,
    skip_transcripts_with_redundant_feature_annotation: bool = True
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


def calculate_optc(
    n_occurrences: int,
    n_counts: int
) -> float:

    return 0.0 if not all((n_occurrences, n_counts)) else float((n_occurrences * 1000) / n_counts)


def calculate_percent_occurrence(
    f_start: int,
    f_end: int,
    n_occurrences: int,
    read_loci: List[Tuple[int, int]]
) -> float:

    n_reads_overlapping_feature = 0

    for read_start, read_end in read_loci:
        if any([
            f_start <= read_start <= f_end,  # Read starts within feature
            f_start <= read_end <= f_end,  # Read ends within feature
            read_start <= f_start <= read_end and read_start <= f_end <= read_end  # Entire feature is within read
        ]):
            n_reads_overlapping_feature += 1

    return 0.0 if not all(
        (n_occurrences, n_reads_overlapping_feature)
    ) else (n_occurrences * 100) / n_reads_overlapping_feature


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

    # Ugly way to make a set from a list of ranges... is there a more pythonic way to do this?
    e_union_set = set()
    for r in e_ranges:
        e_union_set = e_union_set.union(set(r))

    e_j_intersect = e_union_set.intersection(j_range)

    e_f_intersect = e_union_set.intersection(f_range)
    length_of_feature = len(e_f_intersect)

    exonic_overlap = e_j_intersect.intersection(e_f_intersect)
    length_of_overlap = len(exonic_overlap)

    overlap_fraction = 0.0 if not length_of_overlap else length_of_overlap / length_of_feature

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
    read_loci: List[Tuple[int, int]],
    overlap_threshold: float,
    f_start_genomic: int,
    f_end_genomic: int,
    experimental_allow_junction_start_outside_transcript: bool = False,
    experimental_allow_junction_end_outside_transcript: bool = False
) -> OverlappingJunctionsInfo:

    number_occurrences = 0

    overlapping_j_buffer = f""

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

    # Calculate frequency as percent occurrence (splice junctions occur in N% of reads containing part of the feature)
    frequency_occurrences = calculate_percent_occurrence(
        f_start_genomic,
        f_end_genomic,
        number_occurrences,
        read_loci
    )

    return OverlappingJunctionsInfo(
        sample.name,
        number_occurrences,
        frequency_occurrences,
        overlapping_j_buffer
    )


def _convert_loci_to_junction_ids(
    junction_loci: str,
    junction_id_by_locus: Dict[str, str]
) -> str:

    converted_loci = []
    for quantified_locus in junction_loci.split(" "):
        if not quantified_locus:
            continue
        asterisk_index = quantified_locus.index("*")
        quantifier = quantified_locus[:asterisk_index]
        locus = quantified_locus[asterisk_index+2:-1]  # Strip quantifier and square brackets
        converted_loci.append(f"{quantifier}*{junction_id_by_locus[locus]}")

    return " ".join(converted_loci)


def perform_splice_analysis(
    run_name: str,
    features_to_analyze: List[str],
    overlap_threshold: float,
    samples: Dict[str, Sample],
    transcript_library: TranscriptLibrary,
    output_dir: str,
    max_n_features_in_transcript: Union[None, int] = None,
    max_n_processes: int = 1,
    min_bam_parse_time_ns_to_use_multiprocessing: int = 100_000_000,
    min_total_n_junctions_to_use_multiprocessing: int = 10_000,
    mapq_for_unique_mapping: int = DEFAULT_MAPQ_FOR_UNIQUE_MAPPING,
    primary_alignment_only: bool = False,
    include_all_junctions_in_output: bool = True,
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
        print_if_verbose(f"Created worker pool with {n_processes} processes for BAM file processing\n")
    else:
        pool = None

    for feature_substring in features_to_analyze:

        # TODO: Refactor to iterate over transcripts & samples, and perform analysis for all feature substrings with
        #       only one BAM file query. i.e.:
        #       get_splice_junctions_from_sample --> call get_splice_analysis_result_for_sample for each feature
        #       (will require refactoring get_splice_analysis_result_for_sample to accept junctions rather than sample)

        if verbose:
            if max_n_features_in_transcript:
                print(
                    f"Limiting analysis to transcripts containing up to {max_n_features_in_transcript} of the same "
                    f"feature of interest\n"
                )

        print_if_verbose(f"Performing splice analysis for term \"{feature_substring}\"...")

        skipped_exon_match_failure = 0

        with open(os.path.join(output_dir, f"{name_output_file(run_name, feature_substring)}.csv"), "w") as f:

            sample_names_alphabetical = sorted(samples.keys())
            samples_alphabetical = [samples[sample_name] for sample_name in sample_names_alphabetical]

            header = [
                cols.TRANSCRIPT_ID,
                cols.GENE_ID,
                cols.GENE_NAME,
                cols.CHROMOSOME,
                cols.STRAND,
                cols.TRANSCRIPT_START,
                cols.EXON_POSITIONS,
                cols.FEATURE_REGION,
                cols.N_FEATURES_IN_TRANSCRIPT,
                cols.FEATURE_NUMBER,
            ] + ([cols.JUNCTION_DEFINITION, cols.ALL_JUNCTIONS] if include_all_junctions_in_output else []) + [
                cols.OVERLAPPING_JUNCTIONS,
                cols.AVG_NUMBER,
                cols.AVG_FREQUENCY
            ] + sample_names_alphabetical + [f"Percent {sample_name}" for sample_name in sample_names_alphabetical]

            out_csv = csv.writer(f)
            out_csv.writerow(header)

            _info_time_bam_read = False
            _info_time_transcripts = False

            _progress = 0
            for transcripts_by_id in transcript_library.get_transcripts_for_all_genes().values():

                if not transcripts_by_id:
                    continue

                transcripts = list(transcripts_by_id.values())

                if verbose:
                    for _progress_plus_i in [_progress + i for i in range(len(transcripts))]:
                        if _progress_plus_i % 1000 == 0:
                            print(f"Progress: {_progress_plus_i}/{transcript_library.number_of_transcripts}")

                if not any(
                    [feature_substring in transcript.gbseq.get_analysis_features().keys() for transcript in transcripts]
                ):
                    _progress += len(transcripts)
                    continue

                # # START DEBUG BLOCK
                # # if transcripts[0].gene_name not in ("Cacna1d", "Cd74", "Gpr6", "Pkd1"):
                # # if transcripts[0].gene_name not in ("H2-Ab1",):
                # #     _progress += len(transcripts)
                # #     continue
                # if _progress > 1000:
                #     break
                # # print("REGION:", f"{transcript.seqname}:{transcript.start}-{transcript.end}")
                # # END OF DEBUG BLOCK

                chromosome = transcripts[0].seqname
                gene_region_start = min([transcript.start for transcript in transcripts])
                gene_region_end = min([transcript.end for transcript in transcripts])

                if _info_time_bam_read:
                    _bam_read_start_ms = time.time_ns() / 1_000_000

                # Test time taken to parse junctions from first BAM file and decide whether parallelization is necessary
                first_bam_start_ns = time.time_ns()
                first_reads, first_junctions = get_splice_junctions_from_sample(
                    samples_alphabetical[0],
                    chromosome,
                    gene_region_start,
                    gene_region_end,
                    primary_alignment_only,
                    mapq_for_unique_mapping
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
                            primary_alignment_only,
                            mapq_for_unique_mapping
                        ] for sample in samples_alphabetical]
                        remaining_res = pool.starmap(
                            get_splice_junctions_from_sample,
                            args_to_pool
                        )
                        remaining_reads = [res[0] for res in remaining_res]
                        remaining_junctions = [res[1] for res in remaining_res]

                    else:
                        remaining_res = [
                            get_splice_junctions_from_sample(
                                sample,
                                chromosome,
                                gene_region_start,
                                gene_region_end,
                                primary_alignment_only,
                                mapq_for_unique_mapping
                            ) for sample in samples_alphabetical[1:]
                        ]
                        remaining_reads = [res[0] for res in remaining_res]
                        remaining_junctions = [res[1] for res in remaining_res]

                else:

                    remaining_reads, remaining_junctions = [], []

                if _info_time_bam_read:
                    _bam_read_finish_ms = time.time_ns() / 1_000_000
                    print(
                        f"INFO: parsed splice junctions from BAM files for gene {transcripts[0].gene_name} in " +
                        f"{_bam_read_finish_ms - _bam_read_start_ms} ms"
                    )

                all_reads = [first_reads] + remaining_reads
                all_junctions = [first_junctions] + remaining_junctions

                reads_by_sample_name = {
                    sample_names_alphabetical[i]: all_reads[i] for i in range(len(sample_names_alphabetical))
                }
                junctions_by_sample_name = {
                    sample_names_alphabetical[i]: all_junctions[i] for i in range(len(sample_names_alphabetical))
                }

                junction_id_by_locus = {}
                if include_all_junctions_in_output:
                    for junctions in all_junctions:
                        for junction in junctions:
                            locus = f"{junction.start}-{junction.end}"
                            try:
                                # Do nothing if already present as key (rather than iterating over all keys)
                                junction_id_by_locus[locus]
                            except KeyError:
                                junction_id = f"J{1 + len(junction_id_by_locus)}"
                                junction_id_by_locus[locus] = junction_id

                junction_definition = "" if not include_all_junctions_in_output else " ".join(
                    [f"{junction_id}=[{locus}]" for (locus, junction_id) in junction_id_by_locus.items()]
                )

                junctions_loci_by_sample_name = {} if not include_all_junctions_in_output else {
                    sample_name: " ".join(
                        [f"{junc.n_unique}*" + junction_id_by_locus[
                            f"{junc.start}-{junc.end}"
                        ] for junc in junctions if (junc.n_unique > 0)]
                    ) for (sample_name, junctions) in junctions_by_sample_name.items()
                }
                all_junctions_loci = "" if not include_all_junctions_in_output else " | ".join(
                    [f"{name}: {loci}" for (name, loci) in junctions_loci_by_sample_name.items()]
                ).replace("  ", " ")

                for transcript in transcripts:

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

                            feature_region = f"{f_start_genomic}-{f_end_genomic}"
                            exon_positions = " ".join([f"{e.start}-{e.end}" for e in transcript.exons])

                            raw_number_occurrences = {}  # Dict[str, int]
                            frequency_occurrences = {}  # Dict[str, float]

                            if all([
                                pool is not None,
                                sum(
                                    [sum([j.n_unique for j in j_list]) for j_list in junctions_by_sample_name.values()]
                                ) > min_total_n_junctions_to_use_multiprocessing
                            ]):
                                args_to_pool = [[
                                    sample,
                                    junctions_by_sample_name[sample.name],
                                    transcript,
                                    reads_by_sample_name[sample.name],
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
                                    reads_by_sample_name[sample.name],
                                    overlap_threshold,
                                    f_start_genomic,
                                    f_end_genomic
                                ) for sample in samples_alphabetical]

                            for result in results:
                                raw_number_occurrences[result.sample_name] = result.number_occurrences
                                frequency_occurrences[result.sample_name] = result.frequency_occurrences

                            overlapping_junctions_loci = " | ".join(
                                [f"{result.sample_name}: {result.junctions_loci}" for result in results]
                            ) if not include_all_junctions_in_output else " | ".join(
                                [f"""{result.sample_name}: {_convert_loci_to_junction_ids(
                                    result.junctions_loci, junction_id_by_locus
                                )}""" for result in results]
                            )

                            # Write to output (one row per feature per transcript)
                            out_n = [
                                raw_number_occurrences[sample_name] for sample_name in sample_names_alphabetical
                            ]
                            out_freq = [
                                "{:.2f}".format(
                                    frequency_occurrences[sample_name]
                                ) for sample_name in sample_names_alphabetical
                            ]
                            out_csv.writerow(
                                [
                                    transcript.transcript_id,
                                    transcript.gene_id,
                                    transcript.gene_name,
                                    transcript.seqname,
                                    transcript.strand,
                                    transcript.start,
                                    exon_positions,
                                    feature_region,
                                    n_features_in_transcript,
                                    feature_index + 1
                                ] +
                                ([junction_definition, all_junctions_loci] if include_all_junctions_in_output else []) +
                                [overlapping_junctions_loci] +
                                ["{:.2f}".format(mean(out_n))] +
                                ["{:.2f}".format(mean(frequency_occurrences.values()))] +
                                out_n +
                                out_freq
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
            f"...finished ({skipped_exon_match_failure} features skipped either due to failed matching between "
            f"reference and annotation exons, or due to annotated position mapping outside the bounds of the "
            f"corresponding reference exon)\n"
        )
