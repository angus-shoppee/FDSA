
# SPLICE JUNCTION INFORMATION FROM BAM FILES

from typing import List, Tuple, NamedTuple
from enum import Flag, auto
from pysam import AlignmentFile

from experiment import Sample


DEFAULT_MAPQ_FOR_UNIQUE_MAPPING = 255

# CIGAR operations
BAM_CMATCH = 0      # M
BAM_CINS = 1        # I
BAM_CDEL = 2        # D
BAM_CREF_SKIP = 3   # N
BAM_CSOFT_CLIP = 4  # S
BAM_CHARD_CLIP = 5  # H
BAM_CPAD = 6        # P
BAM_CEQUAL = 7      # =
BAM_CDIFF = 8       # X


class SamFlag(Flag):
    READ_PAIRED = auto()
    READ_MAPPED_IN_PROPER_PAIR = auto()
    READ_UNMAPPED = auto()
    MATE_UNMAPPED = auto()
    READ_REVERSE_STRAND = auto()
    MATE_REVERSE_STRAND = auto()
    FIRST_IN_PAIR = auto()
    SECOND_IN_PAIR = auto()
    NOT_PRIMARY_ALIGNMENT = auto()
    READ_FAILS_PLATFORM_QUALITY_CHECKS = auto()
    READ_IS_PCR_OR_OPTICAL_DUPLICATE = auto()
    SUPPLEMENTARY_ALIGNMENT = auto()


_BAM_OP_CONSUMES_REFERENCE = {
    0: True,
    1: False,
    2: True,
    3: True,
    4: False,
    5: False,
    6: False,
    7: True,
    8: True
}


def bam_op_consumes_reference(op: int) -> bool:

    # return True if op in (0, 2, 3, 7, 8) else False
    return _BAM_OP_CONSUMES_REFERENCE[op]


NumberOfJunctions = NamedTuple("NumberOfJunctions", [("unique", int), ("multi", int)])


class SpliceJunction:

    region: str
    chromosome: str
    start: int
    end: int
    n_unique: int
    n_multi: int

    def __init__(
        self,
        region: str,
        n_unique: int,
        n_multi: int
    ) -> None:

        self.region = region
        self.chromosome = region.split(":")[0]
        self.start = int(float(region.split(":")[1].split("-")[0]))
        self.end = int(float(region.split(":")[1].split("-")[1]))

        self.n_unique = n_unique
        self.n_multi = n_multi


def get_splice_junctions_from_sample(
    sample: Sample,
    chromosome: str,
    start: int,
    end: int,
    primary_alignment_only: bool,
    mapq_for_unique_mapping: int,
    offset: int = 0
) -> Tuple[List[Tuple[int, int]], List[SpliceJunction]]:

    sam = AlignmentFile(sample.bam_path, "rb")

    # Fetch reads at the specified region, and extract splice junctions from CIGAR strings

    j_counts = {}  # Dict[str, NumberOfJunctions]

    read_loci = []  # List[Tuple[int, int]]

    for read in sam.fetch(chromosome, start, end):

        if primary_alignment_only:
            if SamFlag.NOT_PRIMARY_ALIGNMENT in SamFlag(read.flag):
                continue

        g = read.reference_start + 1 + offset  # Add 1 as SAM format counts from 1, not 0
        read_start = g

        for op, n in read.cigartuples:  # Iterate over each event/operation in the CIGAR string

            if op == BAM_CREF_SKIP:  # A skip in position on the reference indicates a splice junction
                sj_loc = f"{read.reference_name}:{g}-{g + n - 1}"
                _prev_u, _prev_m = j_counts.get(sj_loc, NumberOfJunctions(0, 0))
                if read.mapping_quality >= mapq_for_unique_mapping:  # Uniquely mapped read
                    j_counts[sj_loc] = NumberOfJunctions(_prev_u + 1, _prev_m)
                else:  # Non-uniquely mapped read
                    j_counts[sj_loc] = NumberOfJunctions(_prev_u, _prev_m + 1)

            if bam_op_consumes_reference(op):  # If this operation "moves along" the reference
                g += n  # Increase genomic position by the given base count

        read_loci.append((read_start, g))

    # Populate junctions

    return (
        read_loci,
        [SpliceJunction(
            region=sj_loc,
            n_unique=num_junctions.unique,
            n_multi=num_junctions.multi
        ) for sj_loc, num_junctions in j_counts.items()]
    )
