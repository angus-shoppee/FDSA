
# SPLICE JUNCTION INFORMATION FROM BAM FILES

from typing import List, NamedTuple
from pysam import AlignmentFile


STAR_DEFAULT_MAPQ_FOR_UNIQUE_MAPPING = 255

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


def bam_op_consumes_reference(op: int) -> bool:

    return True if op in (0, 2, 3, 7, 8) else False


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


def get_splice_junctions_from_sam(
        sam: AlignmentFile,
        chromosome: str,
        start: int,
        end: int,
        offset: int = 0,
) -> List[SpliceJunction]:

    # Fetch reads at the specified region, and extract splice junctions from CIGAR strings

    j_counts = {}

    for read in sam.fetch(chromosome, start, end):

        g = read.reference_start + 1 + offset  # Add 1 as SAM format counts from 1, not 0

        for op, n in read.cigartuples:  # Iterate over each event/operation in the CIGAR string

            if op == BAM_CREF_SKIP:  # A skip in position on the reference indicates a splice junction
                sj_loc = f"{read.reference_name}:{g}-{g + n - 1}"
                _prev_u, _prev_m = j_counts.get(sj_loc, NumberOfJunctions(0, 0))
                if read.mapping_quality >= STAR_DEFAULT_MAPQ_FOR_UNIQUE_MAPPING:  # Uniquely mapped read
                    j_counts[sj_loc] = NumberOfJunctions(_prev_u + 1, _prev_m)
                else:  # Non-uniquely mapped read
                    j_counts[sj_loc] = NumberOfJunctions(_prev_u, _prev_m + 1)

            if bam_op_consumes_reference(op):  # If this operation "moves along" the reference
                g += n  # Increase genomic position by the given base count

    # Populate junctions

    return [SpliceJunction(
        region=sj_loc,
        n_unique=num_junctions.unique,
        n_multi=num_junctions.multi
    ) for sj_loc, num_junctions in j_counts.items()]
