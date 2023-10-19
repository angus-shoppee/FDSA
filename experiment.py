
# EXPERIMENT DESIGN, SAMPLES AND THEIR ASSOCIATED DATA

import os
from pysam import AlignmentFile


class Sample:

    name: str
    suffix: str
    bam_path: str

    def __init__(
        self,
        bam_path: str,
        bam_suffix: str,
    ) -> None:

        self.bam_path = bam_path.replace("\\", "/")
        self.name = os.path.basename(self.bam_path).replace(bam_suffix, "")
        self.suffix = bam_suffix
