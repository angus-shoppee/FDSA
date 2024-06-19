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


# TRANSCRIPT DATA STRUCTURES & ASSOCIATED FUNCTIONS

from dataclasses import dataclass
from typing import List, Dict, Union
from pandas import DataFrame
import os
import pickle

from analysis.annotation import GBSeq, RefseqExon
from analysis.biomart import NameLookup


@dataclass
class ExonRecord:
    exon_number: str  # Starts at 1
    start: int
    end: int
    length: int
    score: str = ""
    sequence: str = ""


def get_exons_from_gtf(
    t_id: str,
    gtf: DataFrame
) -> List[ExonRecord]:

    gtf_transcript_id = list(gtf['transcript_id'])
    gtf_feature = list(gtf['feature'])
    gtf_start = list(gtf['start'])
    gtf_end = list(gtf['end'])
    gtf_score = list(gtf['score'])
    gtf_exon_number = list(gtf['exon_number'])

    exons_feature_indexes = []  # For simplicity
    for i, v in enumerate(gtf_transcript_id):
        if v == t_id and gtf_feature[i] == 'exon':
            exons_feature_indexes.append(i)

    exons = [
        ExonRecord(
            start=gtf_start[i],
            end=gtf_end[i],
            length=(1 + abs(int(gtf_end[i]) - int(gtf_start[i]))),
            score=gtf_score[i],
            exon_number=gtf_exon_number[i]
        ) for i in exons_feature_indexes
    ]
    return exons


def create_exon_map_to_gbseq(
    transcript_exons: List[ExonRecord],
    gbseq_exons: List[RefseqExon]
) -> Dict[str, int]:

    exon_map = {}

    # First, determine the offset needed
    offset = 0
    _n_transcript_exons = len(transcript_exons)
    _n_gbseq_exons = len(gbseq_exons)
    if _n_transcript_exons == 0 or _n_gbseq_exons == 0:
        raise ValueError(
            f"WARNING: Transcript and/or gbseq has zero exons!\n" +
            f"transcript_exons: {[(e.start, e.end) for e in transcript_exons]}\n" +
            f"gbseq_exons: {[(e.start, e.end) for e in gbseq_exons]}"
        )

    if _n_transcript_exons != _n_gbseq_exons:
        _n_find_offset_iters = 1 + abs(_n_transcript_exons - _n_gbseq_exons)
        _with_more_exons = transcript_exons if _n_transcript_exons > _n_gbseq_exons else gbseq_exons
        _with_fewer_exons = transcript_exons if _n_transcript_exons < _n_gbseq_exons else gbseq_exons

        _n_matches_per_round = []
        for _offset in range(_n_find_offset_iters):
            _n_matches = 0
            for i in range(len(_with_fewer_exons)):
                if _with_fewer_exons[i].length == _with_more_exons[i+_offset].length:
                    _n_matches += 1
                _n_matches_per_round.append(_n_matches)

        # # DEBUG BLOCK
        # try:
        #     _n_matches_per_round.index(max(_n_matches_per_round))
        # except ValueError:
        #     print("transcript_exons:", transcript_exons)
        #     print("gbseq_exons:", gbseq_exons)
        #     print("_n_find_offset_iters:", _n_find_offset_iters)
        #     print("_n_matches_per_round:", _n_matches_per_round)
        # # END DEBUG BLOCK

        offset = _n_matches_per_round.index(max(_n_matches_per_round))

    # Then, match individual exons by length and add to map
    # TODO: Optional fallback matching method(s) when lengths are not equal
    # Two cases must be handled separately to avoid wrapping with a negative index

    # Case A: gbseq has more exons
    #         i.e. offset needs to be applied to gbseq_exons
    if _n_transcript_exons < _n_gbseq_exons:
        for i in range(_n_transcript_exons):

            # First, check for overhang
            # Such as in the special case of [a, b, c, d] vs [c, d, e]
            if i + offset > _n_gbseq_exons - 1:
                exon_map[transcript_exons[i].exon_number] = None
                continue

            # Compare (semi-redundantly) exon lengths
            if transcript_exons[i].length == gbseq_exons[i + offset].length:
                exon_map[transcript_exons[i].exon_number] = i + offset
            else:
                exon_map[transcript_exons[i].exon_number] = None

            # The last exon is a special case and can be mapped despite length not matching,
            # as long as the previous exon has also mapped, otherwise don't map as this indicates an edge case
            if i == _n_transcript_exons - 1:
                if _n_transcript_exons > 1:
                    _prev_exon_no = str(int(float(transcript_exons[i].exon_number)) - 1)
                    if exon_map[_prev_exon_no] is not None:
                        exon_map[transcript_exons[i].exon_number] = i + offset

    # Case B: transcript has greater or equal number of exons
    #         i.e. either offset is zero, or needs to be applied to transcript_exons
    else:

        # If offset is applied to transcript_exons, the first element won't map, thus None by default
        if offset > 0:
            exon_map[transcript_exons[0].exon_number] = None

        for i in range(_n_transcript_exons):

            # Skip if (i + offset) is out of bounds
            if i + offset > _n_transcript_exons - 1:
                continue

            # Check for overhang
            elif i > _n_gbseq_exons - 1:
                exon_map[transcript_exons[i + offset].exon_number] = None
                continue

            # Compare lengths
            if transcript_exons[i + offset].length == gbseq_exons[i].length:
                exon_map[transcript_exons[i + offset].exon_number] = i
            else:
                exon_map[transcript_exons[i + offset].exon_number] = None

            # Apply special case for final exon as described above
            if i + offset == _n_transcript_exons - 1:
                if _n_transcript_exons > 1:
                    _prev_exon_no = str(int(float(transcript_exons[i + offset].exon_number)) - 1)
                    if _prev_exon_no in exon_map and exon_map[_prev_exon_no] is not None:
                        exon_map[transcript_exons[i + offset].exon_number] = i

            # TODO: Log transcripts where all keys map to None in exon_map

    return exon_map


class TranscriptRecord:
    """
    Class to store transcript reference features + annotation data
    * Initialize with attributes from a reference GTF
    * Add annotation data by running set_gbseq post-init
    """

    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: str
    strand: str
    frame: str
    transcript_id: str
    gene_name: str
    gene_id: str
    exons: List[ExonRecord]
    gbseq: Union[None, GBSeq]
    exon_map: Dict[str, Union[None, int]]

    def __init__(
        self,
        seqname: str,
        source: str,
        feature: str,
        start: int,
        end: int,
        score: str,
        strand: str,
        frame: str,
        transcript_id: str,
        gene_name: str,
        gene_id: str,
        exons: List[ExonRecord]
    ) -> None:

        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.transcript_id = transcript_id
        self.gene_name = gene_name
        self.gene_id = gene_id

        self.exons = exons

        self.gbseq = None
        self.exon_map = {}  # Exon number (str) mapping to index of refseq exon (int) or None

    def get_length(self) -> int:

        return self.end - self.start if self.end > self.start else self.start - self.end

    def set_gbseq(
        self,
        gbseq: GBSeq
    ) -> None:

        self.gbseq = gbseq

        self.exon_map = create_exon_map_to_gbseq(self.exons, gbseq.exons)

    def genomic_position_from_local(
        self,
        local_position: int,
    ) -> int:

        # TODO: Ensure that positions are being calculated correctly for first/last exons where length doesn't match

        if not self.gbseq:
            raise ValueError("The gbseq attribute must be set first")

        fwd_stranded = True if self.strand == "+" else False

        # Find the exon the local position falls within
        within_gb_exon_index = None
        within_gb_exon_length = None
        distance = local_position
        for i, gbseq_exon in enumerate(self.gbseq.exons):
            if gbseq_exon.start <= local_position <= gbseq_exon.end:
                within_gb_exon_index = i
                within_gb_exon_length = gbseq_exon.length
                break
            distance -= gbseq_exon.length

        inverted_map = {v: k for k, v in self.exon_map.items()}

        _error_exon_info = f"Refseq exons: {[(e.start, e.end) for e in self.exons]}\n" +\
            f"GBSeq exons: {[(e.start, e.end) for e in self.gbseq.exons]}\n" +\
            f"Exon map: {self.exon_map}"

        try:
            within_ref_exon_number = inverted_map[within_gb_exon_index]
            for i, refseq_exon in enumerate(self.exons):
                if int(refseq_exon.exon_number) == int(within_ref_exon_number):
                    delta_length = (refseq_exon.end + 1 - refseq_exon.start) - within_gb_exon_length
                    # Forward stranded: count forward from start of exon
                    # Reverse stranded: count backward from end of exon (due to descending left-to-right formatting)
                    if fwd_stranded:
                        position = refseq_exon.start + distance + delta_length
                    else:
                        position = refseq_exon.end - distance - delta_length
                    if not refseq_exon.start <= position <= refseq_exon.end:
                        raise ValueError(
                            f"Converted local position {local_position} to coordinate {position} which is outside of " +
                            f"the refseq exon (forbidden)\n" + _error_exon_info
                        )
                    return position

        except KeyError:
            raise ValueError(
                f"Failed to convert local position {local_position} - couldn't match exons.\n" + _error_exon_info
            )


class TranscriptLibrary:

    species: str
    number_of_transcripts: int
    _transcripts_by_gene: Dict[str, Dict[str, TranscriptRecord]]

    def __init__(
        self,
        species: str,
        name_lookup: NameLookup,
        ref_gtf: DataFrame,
        numexpr_max_threads: int = 8
    ) -> None:

        # TODO: Replace reference gtf DataFrame with dictionary structure

        self.species = species
        self.number_of_transcripts = 0
        self._transcripts_by_gene = {}

        os.environ["NUMEXPR_MAX_THREADS"] = str(numexpr_max_threads)

        gene_names = name_lookup.get_all_gene_names()

        print("Creating transcript library...")
        _i = 0
        _t = len(gene_names)

        for gene_name in gene_names:

            if _i % 1000 == 0:
                print(f"Progress: ({_i}/{_t})")
            _i += 1

            res = ref_gtf.query(f"gene_name=='{gene_name}'")
            transcripts_slice = res.query(f"feature=='transcript'")

            _t_ids = list(transcripts_slice['transcript_id'])
            transcripts = {}
            # for i, t_id in enumerate(_t_ids):
            for i in range(len(_t_ids)):
                t_id = transcripts_slice.transcript_id.values[i]

                # Only create TranscriptRecord instance if the transcript ID has a matching refseq ID
                # TODO: Log number skipped due to no matching refseq
                try:
                    refseq = name_lookup.convert(t_id, "ensembl", "refseq")
                except KeyError:
                    refseq = ""
                if not refseq:
                    continue

                transcripts[t_id] = TranscriptRecord(
                    seqname=transcripts_slice.seqname.values[i],
                    source=transcripts_slice.source.values[i],
                    feature=transcripts_slice.feature.values[i],
                    start=transcripts_slice.start.values[i],
                    end=transcripts_slice.end.values[i],
                    score=transcripts_slice.score.values[i],
                    strand=transcripts_slice.strand.values[i],
                    frame=transcripts_slice.frame.values[i],
                    transcript_id=transcripts_slice.transcript_id.values[i],
                    gene_name=transcripts_slice.gene_name.values[i],
                    gene_id=transcripts_slice.gene_id.values[i],
                    # exons=get_exons_from_gtf(t_id, transcripts_slice.query(f"transcript_id=='{t_id}'"))
                    exons=get_exons_from_gtf(t_id, res.query(f"transcript_id=='{t_id}'"))
                )
                # _exons = get_exons_from_gtf(t_id, res.query(f"transcript_id=='{t_id}'"))
                # print(f"DEBUG: locus: {transcripts_slice.start.values[i]}-{transcripts_slice.end.values[i]}")
                # print(f"DEBUG: exons: {_exons}")
                self.number_of_transcripts += 1

            self._transcripts_by_gene[gene_name] = transcripts

        print("...done\n")

    def get_all_transcripts(self) -> List[TranscriptRecord]:

        all_transcripts = []
        for d in self._transcripts_by_gene.values():
            all_transcripts += list(d.values())

        return all_transcripts

    def get_transcripts_for_gene(
        self,
        gene_name: str
    ) -> List[TranscriptRecord]:

        return list(self._transcripts_by_gene[gene_name].values())

    def get_transcripts_for_all_genes(self) -> Dict[str, Dict[str, TranscriptRecord]]:

        return self._transcripts_by_gene

    def get_transcript(
        self,
        gene_name: str,
        transcript_id: str
    ) -> TranscriptRecord:

        return self._transcripts_by_gene[gene_name][transcript_id]

    def delete_transcript(
        self,
        gene_name: str,
        transcript_id: str
    ) -> None:

        del self._transcripts_by_gene[gene_name][transcript_id]


def create_and_save_transcript_library(
    species: str,
    output_path: str,
    name_lookup: NameLookup,
    ref_gtf: DataFrame,
    numexpr_max_threads: int = 8,
) -> TranscriptLibrary:

    transcript_library = TranscriptLibrary(
        species,
        name_lookup,
        ref_gtf,
        numexpr_max_threads=numexpr_max_threads
    )

    with open(output_path, "wb") as f:
        pickle.dump(transcript_library, f)

    return transcript_library


def annotate_and_save_transcript_library(
    transcript_library: TranscriptLibrary,
    name_lookup: NameLookup,
    gbseq_by_refseq: Dict[str, GBSeq],
    output_path: str,
) -> TranscriptLibrary:

    """"Modifies transcript_library to add GBSeq objects for transcripts with a refseq_mrna ID;
    all transcripts without a refseq_mrna ID are deleted"""

    # Get all transcripts that have a refseq ID
    print("Applying annotation to transcripts with refseq IDs...")
    _key_errors, _without_refseq, _no_gbseq, _exon_map_failure, _successes = 0, 0, 0, 0, 0
    for transcript in transcript_library.get_all_transcripts():
        # TODO: Remove this check and its counting as this is now handled previously when creating transcript library
        try:
            refseq = name_lookup.convert(transcript.transcript_id, "ensembl", "refseq")
        except KeyError:
            refseq = ""
            _without_refseq -= 1  # Will be re-incremented below; subtracting avoids double counting
            _key_errors += 1
        if refseq:
            try:
                gbseq = gbseq_by_refseq[refseq]
            except KeyError:
                gbseq = None
            if gbseq and gbseq.exons:
                # transcript.set_gbseq(gbseq)
                try:
                    transcript.set_gbseq(gbseq)
                    _successes += 1
                except ValueError:
                    _exon_map_failure += 1
            else:
                transcript_library.delete_transcript(
                    transcript.gene_name,
                    transcript.transcript_id
                )
                _no_gbseq += 1
        else:
            transcript_library.delete_transcript(
                transcript.gene_name,
                transcript.transcript_id
            )
            transcript_library.number_of_transcripts -= 1
            _without_refseq += 1

    print(
        f"... finished with {_successes} successes; removed {_without_refseq} transcripts without refseq ID, " +
        f"{_key_errors} unmappable due to missing or unrecognised ensembl ID, {_no_gbseq} without annotation, "
        f"{_exon_map_failure} failures to map reference genome exons to gbseq exons\n"
    )

    with open(output_path, "wb") as f:
        pickle.dump(transcript_library, f)

    return transcript_library


def load_transcript_library_from_file(path: str):

    with open(path, "rb") as f:

        transcript_library = pickle.load(f)
        # transcript_library = RenameUnpickler(f).load()

    return transcript_library
