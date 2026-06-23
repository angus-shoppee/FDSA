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


# INSPECT TOOLING

from typing import List, Optional
from pathlib import Path
import logging
import gc
import csv
import os

from fdsa.config.parse_config import ProgramRunConfig
from fdsa.analysis.transcript import load_transcript_library_from_file, TranscriptRecord
from fdsa.analysis.core import set_analysis_features


ANNOTATED_TRANSCRIPT_LIBRARY_PICKLE_OBJECT_NAME = "annotated_transcript_library.object"


logger = logging.getLogger(__name__)


# set_analysis_features(
#     run_config.feature_name,
#     annotated_transcript_library,
#     only_use_longest_annotated_transcript=run_config.only_use_longest_annotated_transcript,
#     skip_transcripts_with_redundant_feature_annotation=run_config.skip_transcripts_with_redundant_feature_annotation
# )


def inspect_analysis_features(
        species_specific_data_dir: str,
        output: str,
        run_config: ProgramRunConfig
) -> None:

    inspect_annotations(
        species_specific_data_dir,
        output,
        run_config=run_config
    )


def inspect_annotations(
        species_specific_data_dir: str,
        output: str,
        run_config: Optional[ProgramRunConfig] = None
) -> None:

    simulate_runtime = False if run_config is None else True

    # Possible TODO: Check species against internal config instead of allowing any matching species directory?

    if not os.path.exists(species_specific_data_dir):
        raise ValueError(f"Invalid directory: {species_specific_data_dir}")

    annotated_transcript_library_path = Path(species_specific_data_dir, ANNOTATED_TRANSCRIPT_LIBRARY_PICKLE_OBJECT_NAME)
    if not os.path.exists(annotated_transcript_library_path):
        raise ValueError(f"File \"{ANNOTATED_TRANSCRIPT_LIBRARY_PICKLE_OBJECT_NAME}\" not found in specified directory "
                         f"({species_specific_data_dir})")

    if output.endswith("/") or os.path.isdir(output):
        # Allow either an existing directory, or explicit directory path with trailing slash, and use default filename
        os.makedirs(output, exist_ok=True)
        if simulate_runtime:
            assert run_config is not None
            analysis_features_prefix = f"({run_config.feature_name}) "
        else:
            analysis_features_prefix = ""
        output_path = Path(
            output,
            f"{analysis_features_prefix}{Path(species_specific_data_dir).name}_annotated_transcript_library_export.csv"
        )
    else:
        # If not recognized as a directory, enforce .csv extension to prevent ambiguity
        if not output.endswith(".csv"):
            raise ValueError(f"The specified output path ({output}) was not recognised as a directory and does not end "
                             f"with \".csv\". Please specify either a directory, or a full file path including the "
                             f"final \".csv\" extension.")
        output_path = Path(output)

    logger.info("Loading annotated transcript library...")

    annotated_transcript_library = load_transcript_library_from_file(annotated_transcript_library_path)

    logger.info("... done")

    if simulate_runtime:
        assert run_config is not None
        logger.info(f"Setting analysis features for term: \"{run_config.feature_name}\"...")

        set_analysis_features(
            run_config.feature_name,
            annotated_transcript_library,
            only_use_longest_annotated_transcript=run_config.only_use_longest_annotated_transcript,
            skip_transcripts_with_redundant_feature_annotation=run_config.skip_transcripts_with_redundant_feature_annotation
        )

        logger.info(f"... done")

    logger.info(f"Exporting contents to file ({output_path})...")

    transcripts: List[TranscriptRecord] = annotated_transcript_library.get_all_transcripts()
    del annotated_transcript_library
    gc.collect()

    with open(output_path, "w") as f:

        out = csv.writer(f)

        header = ["gene_name", "gene_id", "transcript_id", "transcript_refseq", "annotations"]
        if simulate_runtime:
            header.append("analysis_features")
        out.writerow(header)

        for t in transcripts:

            if not t.gbseq:
                _e = (f"Unexpectedly encountered transcript without .gbseq attribute while parsing annotated "
                      f"transcript library ({annotated_transcript_library_path}).\nTranscript:\n{t}\n{vars(t)}")
                logger.error(_e)
                raise ValueError(_e)

            row = [t.gene_name, t.gene_id, t.transcript_id, t.gbseq.refseq_id, t.gbseq.serialize_features()]
            if simulate_runtime:
                row.append(t.gbseq.serialize_analysis_features())

            out.writerow(row)

    logger.info(f"... done")



