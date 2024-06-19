# Feature-Directed Splice Analysis
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


# EXPERIMENT DESIGN, SAMPLES AND THEIR ASSOCIATED DATA

import os


class Sample:

    name: str
    bam_ending: str
    bam_path: str

    def __init__(
        self,
        bam_path: str,
        bam_ending: str,
    ) -> None:

        self.bam_path = bam_path.replace("\\", "/")
        self.name = os.path.basename(self.bam_path).replace(bam_ending, "")
        self.bam_ending = bam_ending
