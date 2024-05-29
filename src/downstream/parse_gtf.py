
from dataclasses import dataclass
from typing import Optional, Dict


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
