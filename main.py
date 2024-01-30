
from typing import Any, Union
import gc
import os
import argparse
from gtfparse import read_gtf

from src.config.parse_config import FaseInternalConfig, FaseUserConfig, FaseRunConfig
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

# NEW IN V0.0.5: Run name
RUN_NAME = "Test Run"

# NEW IN V0.0.4: Species from which sequencing data originates
SPECIES = "Mouse"
# SPECIES = "Human"

# NEW IN V0.0.4: Path to single-column csv file / Leave empty to default to whole-genome
# NOTE: Currently unused
GENES_TO_ANALYZE = ""

# REFACTORED IN V0.0.4: Which annotated features do we want to analyze splice events for?
# * values: the strings that will be searched for in the genbank feature's qualifiers
FEATURES_TO_ANALYZE = [
    "transmembrane region"
]

# NEW IN V0.0.4: Indicate whether to only use annotation from the longest annotated transcript per gene
ONLY_USE_LONGEST_ANNOTATED_TRANSCRIPT = True

# NEW IN V0.0.4: Indicate whether to ignore transcripts with duplicate feature annotation
#              (Overridden by ONLY_USE_LONGEST_ANNOTATED_TRANSCRIPT)
SKIP_TRANSCRIPTS_WITH_REDUNDANT_FEATURE_ANNOTATION = True

# NEW IN V0.0.4: Specify whether to only analyze transcripts with a maximum number of n features
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
BAM_ENDING = "_Aligned_Sorted.out.bam"
# BAM_ENDING = "_Aligned.sortedByCoord.out.bam"

# Reference gtf file - make sure these are the same assembly your BAMs were aligned to
REFERENCE_GENOME_GTF_PATH = "/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCm39/gencode_M27_primary/gencode.vM27.primary_assembly.annotation.gtf"
# REFERENCE_GENOME_GTF_PATH = "/Users/aasho2/PHD/Bioinformatics/STAR/genomes/GRCh38/gencode_38_primary/gencode.v38.primary_assembly.annotation.gtf"

# NEW IN V0.0.4: Number of processes to use for reading from multiple bam files simultaneously during splice analysis
# SPLICE_ANALYSIS_MAX_NUMBER_OF_PROCESSES = 24

# NEW IN V0.0.4: Number of threads to use for gtf queries while generating transcript library
# NUMEXPR_MAX_THREADS = 32

# V0.0.5: Replaced NUMEXPR_MAX_THREADS and SPLICE_ANALYSIS_MAX_NUMBER_OF_PROCESSES with shared N_THREADS setting
N_THREADS = 32

# NEW IN V0.0.4: Indicate to perform_splice_analysis whether to only consider primary alignments
#              No longer passed to featureCounts
# TODO: Is this always redundant if only considering unique alignments?
PRIMARY_ALIGNMENT_ONLY = True

FEATURE_JUNCTION_OVERLAP_THRESHOLD = 0.5

EXPERIMENTAL_ALLOW_JUNCTION_START_OUTSIDE_TRANSCRIPT = False
EXPERIMENTAL_ALLOW_JUNCTION_END_OUTSIDE_TRANSCRIPT = False

# Optional: Choose a biomart mirror
BIOMART_MIRROR = 'http://asia.ensembl.org'

# NEW IN V0.0.4
FORCE_REGENERATE_TRANSCRIPT_LIBRARY = False
FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY = False

# NEW IN V0.0.4
FORCE_REGENERATE_WHOLE_GENOME_LOOKUP = False

# NEW IN V0.0.4
FORCE_REGENERATE_GENBANK_FEATURE_XML = False
DO_NOT_RESUME_PARTIAL_DOWNLOAD_GENBANK = False


# ==================================================================================================================== #


PROGRAM_DESCRIPTION = "Feature-directed Analysis of Splice Events (FASE)"

BUILD_NOT_COMPLETED_MESSAGE = ("Annotated transcript library has not yet been built for this species. Please complete "
                               "the build process using \"fase build [species]\"")

INVALID_CONFIG_PATH_MESSAGE = ("The specified path does not represent a valid config file. Please check that the path "
                               "is correct.")

CONFIG_FILE_FORMATTING_ERROR_MESSAGE = ("Error encountered while parsing config file. Please see below for more "
                                        "information:")

USER_CONFIG_NOT_SPECIFIED_MESSAGE = ("A user config file has not been set. Please specify one with "
                                     "\"fase user [path-to-config-file]\"")


# ==================================================================================================================== #

# TODO: Refactor print_if_verbose here and in dependencies to regular if ... -> print
# TODO: Refactor verbose flag to silent flag


def save_user_config_path(
    user_config_path: str,
    save_path_to: str
) -> None:

    if not os.path.exists(user_config_path):
        print(INVALID_CONFIG_PATH_MESSAGE)
        exit()

    try:
        FaseUserConfig(user_config_path)
    except ValueError as e:
        print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
        exit()

    with open(save_path_to, "w") as f:
        f.write(user_config_path)

    print("User config stored.\n")


def build(
    species_specific_data_dir: str,
    species: str,
    genome_gtf_path: str,
    internal_config: FaseInternalConfig,
    verbose: bool = True
):
    def print_if_verbose(*args: Any):
        if verbose:
            print(*args)

    print_if_verbose("Beginning build process...\n")

    # Create data directory for this species if required
    if not os.path.isdir(species_specific_data_dir):
        os.makedirs(species_specific_data_dir)

    # Create name lookup if required, otherwise load from file
    name_lookup_path = os.path.join(species_specific_data_dir, "name_lookup.json")
    if FORCE_REGENERATE_WHOLE_GENOME_LOOKUP or not os.path.exists(name_lookup_path):
        name_lookup = create_and_save_name_lookup(
            internal_config.biomart_name_for_species[SPECIES],
            BIOMART_MIRROR,
            name_lookup_path,
            verbose=verbose
        )
    else:
        name_lookup = load_name_lookup_from_file(name_lookup_path)

    # Create annotated transcript library if required
    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")
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
                numexpr_max_threads=N_THREADS,
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
            os.path.exists(xml_download_progress_file_path),
            # Indicates previously started download has not finished
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
        annotate_and_save_transcript_library(
            transcript_library,
            name_lookup,
            gbseq_by_refseq,
            annotated_transcript_library_path,
            verbose=verbose
        )  # Possible TODO: Remove returned object from this function if no longer used for build-and-run?

        del transcript_library
        del gbseq_by_refseq
        gc.collect()

        print_if_verbose("Annotated transcript library has been created and saved.\n")

    print_if_verbose("Build process complete.\n")


def run(
    species_specific_data_dir: str,
    internal_config: FaseInternalConfig,
    run_config: FaseRunConfig,
    verbose: bool = True
) -> None:

    def print_if_verbose(*args: Any):
        if verbose:
            print(*args)

    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")
    if not os.path.exists(annotated_transcript_library_path):
        print(BUILD_NOT_COMPLETED_MESSAGE)
        return

    # Load annotated transcript library
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

    _ending_length = len(BAM_ENDING)
    bam_file_absolute_paths = [os.path.join(
        BAM_FILES_DIR, p
    ) for p in os.listdir(BAM_FILES_DIR) if p[-_ending_length:] == BAM_ENDING]

    if len(bam_file_absolute_paths) == 0:
        raise ValueError(
            "No BAM files were loaded - check that the supplied BAM directory is valid, and that the specified BAM " +
            f"ending is correct for these files.\nSpecified BAM ending: {BAM_ENDING}\nDetected BAM files: " +
            f"{bam_file_absolute_paths}"
        )

    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, BAM_ENDING)
        samples[sample.name] = sample

    # Run analysis
    perform_splice_analysis(
        RUN_NAME,
        FEATURES_TO_ANALYZE,
        FEATURE_JUNCTION_OVERLAP_THRESHOLD,
        samples,
        annotated_transcript_library,
        OUTPUT_DIR,
        max_n_features_in_transcript=MAX_N_FEATURES_IN_TRANSCRIPT,
        max_n_processes=N_THREADS - 1,
        mapq_for_unique_mapping=internal_config.default_mapq_for_unique_mapping,
        primary_alignment_only=PRIMARY_ALIGNMENT_ONLY,
        include_all_junctions_in_output=INCLUDE_ALL_JUNCTIONS_IN_OUTPUT,
        verbose=verbose
    )


def report() -> None:

    # Generate report
    pass


def main(
    verbose: bool = True
) -> None:

    base_dir = os.path.realpath(os.path.dirname(__file__))
    stored_user_config_path = os.path.join(base_dir, "data", ".user_config_path")

    mode_parser = argparse.ArgumentParser(description=PROGRAM_DESCRIPTION)
    mode_parser.add_argument(
        "mode",
        type=str,
        choices=["user", "build", "run", "report"],
        help="\n".join([
            "user - (Required only once) Sets the user config file containing build information and optional default "
            "parameters.",
            "build - (Required once per species, before first run) Parses transcripts from the provided reference "
            "genome GTF, downloads feature annotation data, and links the two to enable directed splice analysis.",
            "run - Performs splice analysis according to parameters set in the specified run config file. A report "
            "is automatically generated by default, unless otherwise specified.",
            "report - Generates a graphical report from the output of \"fase run\"."
        ])
    )

    mode_arg, subsequent_args = mode_parser.parse_known_args()

    if mode_arg.mode == "user":

        user_parser = argparse.ArgumentParser()
        user_parser.add_argument(
            "user_config_path",
            type=str,
            help="Path to a user config file containing build information and optional default parameters."
        )

        user_args = user_parser.parse_args(subsequent_args)

        save_user_config_path(
            os.path.abspath(user_args.user_config_path),
            stored_user_config_path
        )
        exit()

    else:

        # TODO: if running "fase build", check that species is valid and exit with message if not

        if not os.path.isdir(os.path.join(base_dir, "data")):
            os.mkdir(os.path.join(base_dir, "data"))

        if os.path.exists(stored_user_config_path):
            with open(stored_user_config_path, "r") as f:
                user_config_path = f.read()
            if not os.path.exists(user_config_path):
                print(INVALID_CONFIG_PATH_MESSAGE)
                exit()
        else:
            print(USER_CONFIG_NOT_SPECIFIED_MESSAGE)
            exit()

        internal_config = FaseInternalConfig(
            os.path.join(base_dir, "src", "config", "internal.config")
        )
        try:
            user_config = FaseUserConfig(
                user_config_path
            )
        except ValueError as e:
            print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
            exit()

        if mode_arg.mode == "build":

            build_parser = argparse.ArgumentParser()
            build_parser.add_argument(
                "run_config_path",
                type=str,
                nargs="?",
                default=None,
                help=(
                    "Path to a run config file containing a [RUN] section with the \"species\" and \"genome\" "
                    "parameters set. If not specified, the --species and --genome options required."
                )
            )
            build_parser.add_argument(
                "-s",
                "--species",
                type=str,
                required=False,
                help=(
                    "Species for which transcript annotation data will be downloaded. Must match the species of the "
                    "provided reference genome GTF."
                )
            )
            build_parser.add_argument(
                "-g",
                "--genome",
                type=str,
                required=False,
                help="Reference genome GTF from which transcripts will be parsed."
            )

            build_args = build_parser.parse_args(subsequent_args)

            if build_args.run_config_path is None:
                if not all([
                    build_args.species,
                    build_args.genome
                ]):
                    print("If no run config file is provided, the --species and --genome options must be set.\n")
                    exit()

            if build_args.run_config_path is not None:
                if not os.path.exists(build_args.run_config_path):
                    print(INVALID_CONFIG_PATH_MESSAGE)
                    exit()
                try:
                    run_config = FaseRunConfig(
                        build_args.run_config_path,
                        internal_config,
                        user_config
                    )
                except ValueError as e:
                    print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
                    exit()

                species = build_args.species if build_args.species else run_config.species
                genome = build_args.genome if build_args.genome else run_config.genome

            else:
                # Have previously validated that these are not None if build_args.run_config_path is not set
                species = build_args.species
                genome = build_args.genome

            if species not in internal_config.allowed_species:
                print(f"Invalid species (\"{species}\") specified. Options are: "
                      f"{'/'.join(internal_config.allowed_species)}")
                exit()

            if not os.path.exists(genome):
                print(f"The specified reference genome GTF file was not found. Path: {genome}")

            species_specific_data_dir = os.path.join(base_dir, "data", str(species))

            build(
                species_specific_data_dir,
                species,
                genome,
                internal_config,
                verbose=verbose
            )

        elif mode_arg.mode == "run":

            pass

        exit()
        # --- MOVE BELOW

        species_specific_data_dir = os.path.join(base_dir, "data", SPECIES)

        # Load configs
        try:
            run_config = FaseRunConfig(
                run_config_path,
                internal_config,
                user_config
            )
        except ValueError as e:
            print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
            exit()

        run(
            species_specific_data_dir,
            internal_config,
            verbose=verbose
        )


if __name__ == "__main__":

    main(verbose=True)
