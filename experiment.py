
# EXPERIMENT DESIGN, SAMPLES AND THEIR ASSOCIATED DATA

import os
from pysam import AlignmentFile


class Sample:

    name: str
    suffix: str
    sam: AlignmentFile

    def __init__(
        self,
        bam_path: str,
        bam_suffix: str,
    ) -> None:

        self.name = os.path.basename(bam_path.replace("\\", "/")).replace(bam_suffix, "")
        self.suffix = bam_suffix

        self.sam = AlignmentFile(bam_path, "rb")
