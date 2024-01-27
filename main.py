
from typing import Any
import gc
import os
from gtfparse import read_gtf

from src.config.parse_config import FaseInternalConfig
from src.analysis.transcript import (
    create_and_save_transcript_library, annotate_and_save_transcript_library, load_transcript_library_from_file
)
from src.analysis.entrez import download_and_save_feature_annotation_xml, get_gbseq_from_xml
from src.analysis.biomart import create_and_save_name_lookup, load_name_lookup_from_file
from src.analysis.experiment import Sample
from src.analysis.core import set_analysis_features, perform_splice_analysis

# NOTE: Cd80 does not appear in output but has feature annotation - why?

# ==================================================================================================================== #

# GLOBAL SETTINGS


# Provide user identity to make entrez queries
EMAIL = "angus.shoppee@monash.edu"

# NEW IN V0.4: Species from which sequencing data originates
SPECIES = "Mouse"
# SPECIES = "Human"

# NEW IN V0.4: Path to single-column csv file / Leave empty to default to whole-genome
# NOTE: Currently unused
GENES_TO_ANALYZE = ""

# REFACTORED IN V0.4: Which annotated features do we want to analyze splice events for?
# * values: the strings that will be searched for in the genbank feature's qualifiers
FEATURES_TO_ANALYZE = [
    "transmembrane region"
]

# NEW IN V0.4: Indicate whether to only use annotation from the longest annotated transcript per gene
ONLY_USE_LONGEST_ANNOTATED_TRANSCRIPT = True

# NEW IN V0.4: Indicate whether to ignore transcripts with duplicate feature annotation
#              (Overridden by ONLY_USE_LONGEST_ANNOTATED_TRANSCRIPT)
SKIP_TRANSCRIPTS_WITH_REDUNDANT_FEATURE_ANNOTATION = True

# NEW IN V0.4: Specify whether to only analyze transcripts with a maximum number of n features
#              (set to None or 0 to analyze all transcripts regardless of feature number)
MAX_N_FEATURES_IN_TRANSCRIPT = 1

# Where to save output files
OUTPUT_DIR = "/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4"
# OUTPUT_DIR = "/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4 Mutu2020"
# OUTPUT_DIR = "/Users/aasho2/Projects/FASE_V1/OUTPUT/V0_4 HumanBloodDC"

# NEW IN V0.4: Indicate whether to output all splice junctions for each sample alongside overlapping junctions
#              (Required for plotting splice graphs, but will increase output file size)
#              NOTE: If number of unique splice junctions is very high, may exceed the Excel cell character limit
INCLUDE_ALL_JUNCTIONS_IN_OUTPUT = True

# Sorted, indexed BAM files to be analyzed
# BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/hons_PD1KO/sorted"
BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/MuTu_dabraf_2020/sorted"
# BAM_FILES_DIR = "/Users/aasho2/PHD/Bioinformatics/STAR/runs/PRJNA287649_human_blood_dcs/run_231023/star_output"

# Everything that comes after sample identifier in the bam paths, including file extension
BAM_SUFFIX = "_Aligned_Sorted.out.bam"
# BAM_SUFFIX = "_Aligned.sortedByCoord.out.bam"

# Reference gtf file - make sure these are the same assembly your BAMs were aligned to
REFERENCE_GENOME_GTF_PATH = "/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCm39/gencode_M27_primary/gencode.vM27.primary_assembly.annotation.gtf"
# REFERENCE_GENOME_GTF_PATH = "/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCh38/gencode_38_primary/gencode.v38.primary_assembly.annotation.gtf"

# NEW IN V0.4: Number of processes to use for reading from multiple bam files simultaneously during splice analysis
SPLICE_ANALYSIS_MAX_NUMBER_OF_PROCESSES = 24

# NEW IN V0.4: Number of threads to use for gtf queries while generating transcript library
NUMEXPR_MAX_THREADS = 32

# NEW IN V0.4: Indicate to perform_splice_analysis whether to only consider primary alignments
#              No longer passed to featureCounts
# TODO: Is this always redundant if only considering unique alignments?
PRIMARY_ALIGNMENT_ONLY = True

FEATURE_JUNCTION_OVERLAP_THRESHOLD = 0.5

EXPERIMENTAL_ALLOW_JUNCTION_START_OUTSIDE_TRANSCRIPT = False
EXPERIMENTAL_ALLOW_JUNCTION_END_OUTSIDE_TRANSCRIPT = False

# Optional: Choose a biomart mirror
BIOMART_MIRROR = 'http://asia.ensembl.org'

# NEW IN V0.4
FORCE_REGENERATE_TRANSCRIPT_LIBRARY = False
FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY = False

# NEW IN V0.4
FORCE_REGENERATE_WHOLE_GENOME_LOOKUP = False

# NEW IN V0.4
FORCE_REGENERATE_GENBANK_FEATURE_XML = False
DO_NOT_RESUME_PARTIAL_DOWNLOAD_GENBANK = False


# ==================================================================================================================== #

def main(verbose: bool = False) -> None:

    # TODO: Refactor print_if_verbose here and in dependencies to print / print_if_not_silent

    def print_if_verbose(*args: Any):
        if verbose:
            print(*args)

    print_if_verbose("Initializing...")

    base_dir = os.path.realpath(os.path.dirname(__file__))

    # Load config
    fase_analysis_config = FaseInternalConfig(
        os.path.join(base_dir, "src", "config", "internal.config")
    )

    # Validate user settings (not exhaustive)
    # TODO: Complete
    if SPECIES not in fase_analysis_config.allowed_species:
        raise ValueError(
            f"Provided species '{SPECIES}' is not supported. Currently supported species: {fase_analysis_config.allowed_species}"
        )

    # Create data directory for this species if required
    species_specific_data_dir = os.path.join(base_dir, "data", SPECIES)
    if not os.path.isdir(species_specific_data_dir):
        os.makedirs(species_specific_data_dir)

    name_lookup_path = os.path.join(species_specific_data_dir, "name_lookup.json")

    # Create name lookup if required
    if FORCE_REGENERATE_WHOLE_GENOME_LOOKUP or not os.path.exists(name_lookup_path):
        name_lookup = create_and_save_name_lookup(
            fase_analysis_config.biomart_name_for_species[SPECIES],
            BIOMART_MIRROR,
            name_lookup_path,
            verbose=verbose
        )

    # Load name lookup
    else:
        name_lookup = load_name_lookup_from_file(name_lookup_path)

    print_if_verbose("... initialization finished\n")

    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")

    # Create annotated transcript library if required
    if any([
        FORCE_REGENERATE_TRANSCRIPT_LIBRARY,
        FORCE_REGENERATE_GENBANK_FEATURE_XML,
        FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY,
        not os.path.exists(annotated_transcript_library_path)
    ]):

        if not os.path.exists(annotated_transcript_library_path):
            print_if_verbose("Annotated transcript library has not been built and will be initialized this run.\n")
        if FORCE_REGENERATE_TRANSCRIPT_LIBRARY:
            print_if_verbose("The FORCE_REGENERATE_TRANSCRIPT_LIBRARY flag has been enabled.\n")
        if FORCE_REGENERATE_GENBANK_FEATURE_XML:
            print_if_verbose("The FORCE_REGENERATE_GENBANK_FEATURE_XML flag has been enabled.\n")
        if FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY:
            print_if_verbose("The FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY flag has been enabled.\n")

        transcript_library_path = os.path.join(species_specific_data_dir, "transcript_library.object")

        # Create transcript library if required
        if FORCE_REGENERATE_TRANSCRIPT_LIBRARY or not os.path.exists(transcript_library_path):
            print_if_verbose("Importing reference GTF...")
            ref_gtf = read_gtf(REFERENCE_GENOME_GTF_PATH)
            print_if_verbose("...done\n")

            transcript_library = create_and_save_transcript_library(
                SPECIES,
                transcript_library_path,
                name_lookup,
                ref_gtf,
                numexpr_max_threads=NUMEXPR_MAX_THREADS,
                verbose=verbose
            )
            del ref_gtf
            gc.collect()

        # Load transcript library
        else:
            transcript_library = load_transcript_library_from_file(transcript_library_path)

        genbank_feature_xml_path = os.path.join(species_specific_data_dir, "genbank_feature_annotation.xml")
        xml_download_progress_file_path = os.path.join(species_specific_data_dir, ".download_progress.json")

        # Create genbank feature XML if required
        if any([
            FORCE_REGENERATE_GENBANK_FEATURE_XML,
            os.path.exists(xml_download_progress_file_path),  # Indicates previously started download has not finished
            not os.path.exists(genbank_feature_xml_path)
        ]):
            print_if_verbose("Downloading genbank feature annotation data...")

            transcripts_with_refseq = []
            for t in transcript_library.get_all_transcripts():
                try:
                    refseq = name_lookup.convert(t.transcript_id, "ensembl", "refseq")
                    if refseq:
                        # transcripts_with_refseq.append(t.transcript_id)
                        transcripts_with_refseq.append(refseq)
                except KeyError:
                    pass

            download_and_save_feature_annotation_xml(
                list(transcripts_with_refseq),
                genbank_feature_xml_path,
                EMAIL,
                force_restart_download=DO_NOT_RESUME_PARTIAL_DOWNLOAD_GENBANK,
                verbose=verbose
            )

            del transcripts_with_refseq
            gc.collect()

            print_if_verbose("...done\n")

        # Load genbank feature XML
        gbseq_by_refseq = get_gbseq_from_xml(
            genbank_feature_xml_path,
            verbose=verbose
        )

        # Apply feature annotation to transcript library
        annotated_transcript_library = annotate_and_save_transcript_library(
            transcript_library,
            name_lookup,
            gbseq_by_refseq,
            annotated_transcript_library_path,
            verbose=verbose
        )

        del transcript_library
        del gbseq_by_refseq
        gc.collect()

        print_if_verbose("Annotated transcript library has been created and saved.\n")

    # Load annotated transcript library
    else:
        print_if_verbose("Loading transcript library...")
        annotated_transcript_library = load_transcript_library_from_file(annotated_transcript_library_path)
        print_if_verbose("...done\n")

    print_if_verbose(f"Enabling screening for features containing term(s): \"{', '.join(FEATURES_TO_ANALYZE)}\"...")

    set_analysis_features(
        FEATURES_TO_ANALYZE,
        annotated_transcript_library,
        only_use_longest_annotated_transcript=ONLY_USE_LONGEST_ANNOTATED_TRANSCRIPT,
        skip_transcripts_with_redundant_feature_annotation=SKIP_TRANSCRIPTS_WITH_REDUNDANT_FEATURE_ANNOTATION
    )

    print_if_verbose("... done\n")

    if not os.path.isdir(OUTPUT_DIR):
        if os.path.exists(OUTPUT_DIR):
            raise ValueError(f"The specified output directory ({OUTPUT_DIR}) already exists as a file.")
        os.makedirs(OUTPUT_DIR)

    _suffix_length = len(BAM_SUFFIX)
    bam_file_absolute_paths = [os.path.join(
        BAM_FILES_DIR, p
    ) for p in os.listdir(BAM_FILES_DIR) if p[-_suffix_length:] == BAM_SUFFIX]

    if len(bam_file_absolute_paths) == 0:
        raise ValueError(
            "No BAM files were loaded - check that the supplied BAM directory is valid, and that the specified BAM " +
            f"suffix is correct for these files.\nSpecified BAM suffix: {BAM_SUFFIX}\nDetected BAM files: " +
            f"{bam_file_absolute_paths}"
        )

    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, BAM_SUFFIX)
        samples[sample.name] = sample

    # Run analysis
    perform_splice_analysis(
        FEATURES_TO_ANALYZE,
        FEATURE_JUNCTION_OVERLAP_THRESHOLD,
        samples,
        annotated_transcript_library,
        OUTPUT_DIR,
        max_n_features_in_transcript=MAX_N_FEATURES_IN_TRANSCRIPT,
        max_n_processes=SPLICE_ANALYSIS_MAX_NUMBER_OF_PROCESSES,
        primary_alignment_only=PRIMARY_ALIGNMENT_ONLY,
        include_all_junctions_in_output=INCLUDE_ALL_JUNCTIONS_IN_OUTPUT,
        verbose=verbose
    )


if __name__ == "__main__":

    main(verbose=True)
