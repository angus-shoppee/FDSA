
from dataclasses import dataclass, asdict
from typing import Optional, Tuple, Dict


STRINGTIE_SOURCE_NAME = "StringTie"

SEQNAME_INDEX = 0
SOURCE_INDEX = 1
FEATURE_INDEX = 2
START_INDEX = 3
END_INDEX = 4
SCORE_INDEX = 5
STRAND_INDEX = 6
FRAME_INDEX = 7
ATTRIBUTES_INDEX = 8


@dataclass
class GtfRecord:
    seqname: str
    source: str
    feature: str
    start: int
    end: int
    score: Optional[float]
    strand: Optional[str]
    frame: Optional[int]
    attributes: Dict[str, str]


@dataclass
class GtfTranscript(GtfRecord):
    exons: Dict[str, Tuple[int, int]]

    def __init__(
        self,
        gtf_record: GtfRecord
    ) -> None:

        super().__init__(**asdict(gtf_record))

        self.exons = {}

    def set_exon(
        self,
        exon_number: str,
        exon_start: int,
        exon_end: int
    ) -> None:

        self.exons[exon_number] = (exon_start, exon_end)


def parse_gtf_record(line: str) -> GtfRecord:

    fields = line.split("\t")

    return GtfRecord(
        seqname=fields[SEQNAME_INDEX],
        source=fields[SOURCE_INDEX],
        feature=fields[FEATURE_INDEX],
        start=int(fields[START_INDEX]),
        end=int(fields[END_INDEX]),
        score=None if fields[SCORE_INDEX] == "." else float(fields[SCORE_INDEX]),
        strand=None if fields[STRAND_INDEX] == "." else fields[STRAND_INDEX],
        frame=None if fields[FRAME_INDEX] == "." else int(fields[FRAME_INDEX]),
        attributes={
            pair[0]: pair[1].replace("\"", "")  # Do not include quotation mark characters in tag values
            for pair in [
                s.strip().split(" ")
                for s in fields[ATTRIBUTES_INDEX].split(";")
            ]
            if len(pair) == 2
        }
    )


def get_stringtie_transcripts_from_gtf(gtf_path: str) -> Dict[str, Dict[str, GtfTranscript]]:

    # Transcripts are stored against transcript IDs, collectively stored against gene IDs

    results: Dict[str, Dict[str, GtfTranscript]] = {}

    with open(gtf_path, "r") as gtf_file:

        # Parse all transcripts before parsing any exons, in case an exon record precedes its corresponding transcript
        for line in gtf_file:

            if line[0] == "#":
                continue

            # For a small speed improvement, only parse lines fully if they are StringTie transcripts
            _line_split = line.split("\t")
            if _line_split[SOURCE_INDEX] == STRINGTIE_SOURCE_NAME and _line_split[FEATURE_INDEX] == "transcript":

                transcript_record = parse_gtf_record(line)

                gene_id = transcript_record.attributes["gene_id"]
                transcript_id = transcript_record.attributes["transcript_id"]

                if not results.get(gene_id):
                    results[gene_id] = {}

                results[gene_id][transcript_id] = GtfTranscript(transcript_record)

        # Parse exons and populate corresponding transcripts' exon info
        gtf_file.seek(0)
        for line in gtf_file:

            if line[0] == "#":
                continue

            _line_split = line.split("\t")
            if _line_split[SOURCE_INDEX] == STRINGTIE_SOURCE_NAME and _line_split[FEATURE_INDEX] == "exon":

                exon_record = parse_gtf_record(line)

                gene_id = exon_record.attributes["gene_id"]
                transcript_id = exon_record.attributes["transcript_id"]

                results[gene_id][transcript_id].set_exon(
                    exon_record.attributes["exon_number"],
                    exon_record.start,
                    exon_record.end
                )

    return results
