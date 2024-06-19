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


from typing import Union
import gc
import os
import argparse
from sys import argv
# TODO: Remove gtfparse dependency
from gtfparse import read_gtf

from config.parse_args import (
    get_mode_parser, get_build_parser, get_run_parser, get_report_parser, get_filter_parser, get_quant_parser
)
from config.parse_config import ProgramInternalConfig, ProgramUserConfig, ProgramRunConfig
from analysis.transcript import (
    create_and_save_transcript_library, annotate_and_save_transcript_library, load_transcript_library_from_file
)
from analysis.entrez import download_and_save_feature_annotation_xml, get_gbseq_from_xml
from analysis.biomart import create_and_save_name_lookup, load_name_lookup_from_file
from analysis.experiment import Sample
from analysis.core import set_analysis_features, perform_splice_analysis, name_output
from reporting.report import create_report
from downstream.filter import generate_filtered_bam_files, FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME
from downstream.quant import quantify_isoforms


# ==================================================================================================================== #


PROGRAM_DESCRIPTION = "Feature-Directed Splice Analysis"

FDSA_BUILD_COMMAND_USAGE = ("fdsa build RUN_CONFIG_PATH\n(or)\n"
                            "fdsa build --species SPECIES_NAME --genome REFERENCE_GENOME_GTF_PATH")

FDSA_RUN_COMMAND_USAGE = "fdsa run [--report] [--no-report] RUN_CONFIG_PATH"

FDSA_REPORT_COMMAND_USAGE = "fdsa report RUN_CONFIG_PATH"

FDSA_FILTER_COMMAND_USAGE = "fdsa filter RUN_CONFIG_PATH"

FDSA_QUANT_COMMAND_USAGE = ("fdsa quant RUN_CONFIG_PATH [-g REFERENCE_GENOME_PATH] [-i FILTERED_BAM_DIRECTORY] "
                            "[-o OUTPUT_DIRECTORY]\n(or)\n"
                            "fdsa quant -g REFERENCE_GENOME_PATH -i FILTERED_BAM_DIRECTORY -o OUTPUT_DIRECTORY")

BUILD_NOT_COMPLETED_MESSAGE = ("Annotated transcript library has not yet been built for this species. Please complete "
                               "the build process using \"fdsa build [species]\"")

INVALID_CONFIG_PATH_MESSAGE = ("The specified path does not represent a valid config file. Please check that the path "
                               "is correct.")

CONFIG_FILE_FORMATTING_ERROR_MESSAGE = ("Error encountered while parsing config file. Please see below for more "
                                        "information:")

USER_CONFIG_NOT_SPECIFIED_MESSAGE = ("A user config file has not been set. No user defaults will be loaded, and the "
                                     "\"fdsa build\" command will require an email_for_apis address to be supplied at "
                                     "run time. To set a user config file, run: \"fdsa user [path-to-config-file]\"")

NO_EMAIL_ADDRESS_MESSAGE = ("An email address must be provided in order to make GenBank API requests. Please "
                            "either set a user config with a valid \"email\" parameter under the [BUILD] "
                            "section, or otherwise provide an email_for_apis address using the \"--email\" "
                            "option.")


# ==================================================================================================================== #

# TODO: Refactor print to use logging, and add optional silent mode for progress indicators etc.

# NOTE: Cd80 does not appear in output but has feature annotation - why?


def _confirm_build_overwrite() -> bool:

    user_choice = input("This action will overwrite existing build files. Proceed? "
                        "(enter y/yes to proceed, any other key to exit). ")

    return True if user_choice.lower() in ("y", "yes") else False


def _user(
    user_config_path: str,
    save_path_to: str
) -> None:

    if not os.path.exists(user_config_path):
        print(INVALID_CONFIG_PATH_MESSAGE)
        exit()

    try:
        ProgramUserConfig(user_config_path)
    except ValueError as e:
        print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
        exit()

    with open(save_path_to, "w") as f:
        f.write(user_config_path)

    print("User config stored.\n")


def _build(
    species_specific_data_dir: str,
    species: str,
    genome_gtf_path: str,
    email: str,
    n_threads: int,
    internal_config: ProgramInternalConfig,
    user_config: ProgramUserConfig,
    force_redo_annotate_transcript_library: bool = False,
    force_regenerate_whole_genome_lookup: bool = False,
    force_regenerate_transcript_library: bool = False,
    force_regenerate_genbank_features: bool = False,
    restart_partial_download: bool = False
):

    if any([
        force_redo_annotate_transcript_library,
        force_regenerate_whole_genome_lookup,
        force_regenerate_transcript_library,
        force_regenerate_genbank_features
    ]):
        if not _confirm_build_overwrite():
            exit()

    print("Beginning build process...\n")

    # Create data directory for this species if required
    if not os.path.isdir(species_specific_data_dir):
        os.makedirs(species_specific_data_dir)

    # Create name lookup if required, otherwise load from file
    name_lookup_path = os.path.join(species_specific_data_dir, "name_lookup.json")
    if force_regenerate_whole_genome_lookup or not os.path.exists(name_lookup_path):
        name_lookup = create_and_save_name_lookup(
            internal_config.biomart_name_for_species[species],
            (user_config.user_default_biomart_mirror if user_config.user_default_biomart_mirror is not None
             else internal_config.default_biomart_mirror),
            name_lookup_path
        )
    else:
        name_lookup = load_name_lookup_from_file(name_lookup_path)

    # Create annotated transcript library if required
    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")
    if any([
        force_regenerate_transcript_library,
        force_regenerate_genbank_features,
        force_redo_annotate_transcript_library,
        not os.path.exists(annotated_transcript_library_path)
    ]):

        if not os.path.exists(annotated_transcript_library_path):
            print("Annotated transcript library has not been built and will be initialized this run.\n")
        if force_regenerate_transcript_library:
            print("The FORCE_REGENERATE_TRANSCRIPT_LIBRARY flag has been enabled.\n")
        if force_regenerate_genbank_features:
            print("The FORCE_REGENERATE_GENBANK_FEATURE_XML flag has been enabled.\n")
        if force_redo_annotate_transcript_library:
            print("The FORCE_REDO_ANNOTATE_TRANSCRIPT_LIBRARY flag has been enabled.\n")

        transcript_library_path = os.path.join(species_specific_data_dir, "transcript_library.object")

        # Create transcript library if required
        if force_regenerate_transcript_library or not os.path.exists(transcript_library_path):
            print("Importing reference GTF...")
            # TODO: Use existing internal implementation of GTF parsing and remove gtfparse dependency
            ref_gtf = read_gtf(genome_gtf_path)
            print("...done\n")

            transcript_library = create_and_save_transcript_library(
                species,
                transcript_library_path,
                name_lookup,
                ref_gtf,
                numexpr_max_threads=n_threads
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
            force_regenerate_genbank_features,
            os.path.exists(xml_download_progress_file_path),
            # Indicates previously started download has not finished
            not os.path.exists(genbank_feature_xml_path)
        ]):
            print("Downloading genbank feature annotation data...")

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
                email,
                force_restart_download=restart_partial_download
            )

            del transcripts_with_refseq
            gc.collect()

            print("...done\n")

        # Load genbank feature XML
        gbseq_by_refseq = get_gbseq_from_xml(
            genbank_feature_xml_path
        )

        # Apply feature annotation to transcript library
        annotate_and_save_transcript_library(
            transcript_library,
            name_lookup,
            gbseq_by_refseq,
            annotated_transcript_library_path
        )  # Possible TODO: Remove returned object from this function if no longer used for build-and-run?

        del transcript_library
        del gbseq_by_refseq
        gc.collect()

        print("Annotated transcript library has been created and saved.\n")

    print("Build process complete.\n")


def _run(
    species_specific_data_dir: str,
    internal_config: ProgramInternalConfig,
    run_config: ProgramRunConfig
) -> None:

    annotated_transcript_library_path = os.path.join(species_specific_data_dir, "annotated_transcript_library.object")
    if not os.path.exists(annotated_transcript_library_path):
        print(BUILD_NOT_COMPLETED_MESSAGE)
        return

    # Load annotated transcript library
    print(f"Loading {run_config.species} transcript library...")
    annotated_transcript_library = load_transcript_library_from_file(annotated_transcript_library_path)
    print("...done\n")

    print(f"Enabling screening for features containing term: \"{run_config.feature_name}\"...")

    set_analysis_features(
        run_config.feature_name,
        annotated_transcript_library,
        only_use_longest_annotated_transcript=run_config.only_use_longest_annotated_transcript,
        skip_transcripts_with_redundant_feature_annotation=run_config.skip_transcripts_with_redundant_feature_annotation
    )

    print("... done\n")

    if not os.path.isdir(run_config.output_path):
        if os.path.exists(run_config.output_path):
            raise ValueError(f"The specified output directory ({run_config.output_path}) already exists as a file.")
        os.makedirs(run_config.output_path)

    _ending_length = len(run_config.bam_ending)
    bam_file_absolute_paths = [os.path.join(
        run_config.input_path, p
    ) for p in os.listdir(run_config.input_path) if p[-_ending_length:] == run_config.bam_ending]

    if len(bam_file_absolute_paths) == 0:
        raise ValueError(
            "No BAM files were loaded - check that the supplied BAM directory is valid, and that the specified BAM " +
            f"ending is correct for these files.\nSpecified BAM ending: {run_config.bam_ending}\nDetected BAM files: " +
            f"{bam_file_absolute_paths}"
        )

    samples = {}
    for bam_path in bam_file_absolute_paths:
        sample = Sample(bam_path, run_config.bam_ending)
        samples[sample.name] = sample

    # Run analysis
    perform_splice_analysis(
        run_config.run_name,
        run_config.feature_name,
        run_config.feature_junction_overlap_threshold,
        run_config.genes,
        samples,
        annotated_transcript_library,
        run_config.output_path,
        max_n_features_in_transcript=run_config.max_n_features_in_transcript,
        max_n_processes=run_config.n_threads - 1,
        mapq_for_unique_mapping=internal_config.default_mapq_for_unique_mapping,
        primary_alignment_only=run_config.primary_alignment_only,
        include_all_junctions_in_output=run_config.include_all_junctions_in_output
    )


def _report(
    run_config: ProgramRunConfig
) -> None:

    create_report(run_config)


def _filter(
    run_config: ProgramRunConfig
) -> None:

    generate_filtered_bam_files(run_config)


def _quant(
    run_config: Union[ProgramRunConfig, None],
    user_config: Union[ProgramUserConfig, None] = None,
    stringtie_executable_path: Union[str, None] = None,
    prep_de_script_path: Union[str, None] = None,
    gtf_path: Union[str, None] = None,
    filtered_bam_dir: Union[str, None] = None,
    output_dir: Union[str, None] = None,
    read_length: Union[int, None] = None,
    fdsa_results_path: Union[str, None] = None,
    disable_assign_reference_gene: Union[bool, None] = None,
    threads: Union[int, None] = None  # Allow None to indicate that value should be loaded from run config
) -> None:

    # Can be run either stand-alone (all kwargs required), or with a run config, allowing overrides by defining kwargs

    # Ensure all parameters have been set if both run config and user config are omitted
    if run_config is None:

        if any([gtf_path is None, filtered_bam_dir is None, output_dir is None]):
            raise ValueError("If a run config is not supplied, then all of the gtf_path, filtered_bam_dir, and "
                             "output_dir parameters must be set.")

        if user_config is None and any([
            stringtie_executable_path is None,
            prep_de_script_path is None
        ]):

            raise ValueError(f"Paths to the stringtie executable and prepDE.py3 script are required. These "
                             f"must either be set using command line arguments, or defined in either the "
                             f"[QUANT] section of run config or the [DEFAULT QUANT] section of user config. "
                             f"No user config was found - consider creating one (see fase user --help for "
                             f"more information).")

    if run_config is None:

        _stringtie_executable_path = stringtie_executable_path if stringtie_executable_path is not None \
            else user_config.user_default_quant_stringtie_executable_path
        _prep_de_script_path = prep_de_script_path if prep_de_script_path is not None \
            else user_config.user_default_quant_prep_de_script_path
        _gtf_path = gtf_path
        _filtered_bam_dir = filtered_bam_dir
        _output_dir = output_dir
        _read_length = read_length
        _fdsa_results_path = fdsa_results_path
        _assign_reference_gene = False if disable_assign_reference_gene else True
        _threads = threads

    else:

        inferred_filtered_bam_dir = os.path.abspath(
            os.path.join(
                run_config.output_path,
                FILTERED_BAM_FILES_DEFAULT_DIRECTORY_NAME,
                f"{name_output(run_config.run_name, run_config.feature_name)}"
            )
        )

        inferred_fdsa_results_path = os.path.join(
            run_config.output_path,
            f"{name_output(run_config.run_name, run_config.feature_name)}.csv"
        )

        _stringtie_executable_path = stringtie_executable_path if stringtie_executable_path is not None \
            else run_config.quant_stringtie_executable_path
        _prep_de_script_path = prep_de_script_path if prep_de_script_path is not None \
            else run_config.quant_prep_de_script_path
        _gtf_path = gtf_path if gtf_path is not None else run_config.genome
        _filtered_bam_dir = filtered_bam_dir if filtered_bam_dir is not None else inferred_filtered_bam_dir
        _output_dir = output_dir if output_dir is not None else run_config.output_path
        _read_length = read_length if read_length is not None else run_config.quant_read_length
        _fdsa_results_path = fdsa_results_path if fdsa_results_path is not None else (
            inferred_fdsa_results_path if os.path.isfile(inferred_fdsa_results_path) else None
        )
        _assign_reference_gene = False if disable_assign_reference_gene else run_config.quant_assign_reference_gene
        _threads = threads if threads is not None else run_config.n_threads

    # Check that _stringtie_executable_path and _prep_de_script_path have been defined in one of the three possible ways
    if any([
        _stringtie_executable_path is None,
        _prep_de_script_path is None
    ]):
        raise ValueError(f"The parameter(s) stringtieExecutablePath and/or prepDeScriptPath have not been set "
                         f"either in the [QUANT] section of run config, or the [DEFAULT QUANT] section of user "
                         f"config. If neither are defined, then these parameters must be set using command line "
                         f"arguments (see fdsa quant --help for more information).")

    quantify_isoforms(
        _stringtie_executable_path,
        _prep_de_script_path,
        _gtf_path,
        _filtered_bam_dir,
        _output_dir,
        read_length=_read_length,
        fdsa_results_path=_fdsa_results_path,
        assign_reference_gene=_assign_reference_gene,
        threads=_threads
    )


def main() -> None:

    # TODO: Auto-generated usage at top line of arg parser help text is incorrect for modes, need override

    base_dir = os.path.realpath(os.path.dirname(__file__))
    stored_user_config_path = os.path.join(base_dir, "data", ".user_config_path")

    mode_parser = get_mode_parser(PROGRAM_DESCRIPTION)

    if len(argv) > 1 and argv[1] in ("-h", "--help"):
        mode_parser.print_help()
        exit()

    mode_arg, subsequent_args = mode_parser.parse_known_args()

    if mode_arg.mode == "user":

        user_parser = argparse.ArgumentParser()
        user_parser.add_argument(
            "user_config_path",
            type=str,
            help="Path to a user config file containing build information and optional default parameters."
        )

        user_args = user_parser.parse_args(subsequent_args)

        _user(
            os.path.abspath(user_args.user_config_path),
            stored_user_config_path
        )
        exit()

    else:

        if not os.path.isdir(os.path.join(base_dir, "data")):
            os.mkdir(os.path.join(base_dir, "data"))

        if os.path.exists(stored_user_config_path):
            with open(stored_user_config_path, "r") as f:
                user_config_path = f.read()
            if not os.path.exists(user_config_path):
                print(INVALID_CONFIG_PATH_MESSAGE)
                exit()
        else:
            print("NOTE:", USER_CONFIG_NOT_SPECIFIED_MESSAGE, "\n")
            user_config_path = None

        internal_config = ProgramInternalConfig(
            os.path.join(base_dir, "src", "config", "internal.config")
        )
        if user_config_path is not None:
            try:
                user_config = ProgramUserConfig(
                    user_config_path
                )
            except ValueError as e:
                print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
                exit()
        else:
            user_config = None

        if mode_arg.mode == "build":

            build_parser = get_build_parser()

            if any([arg in subsequent_args for arg in ("-h", "--help")]):
                build_parser.print_help()
                exit()

            build_args = build_parser.parse_args(subsequent_args)

            if build_args.run_config_path is None:
                if not all([
                    build_args.species,
                    build_args.genome,
                    build_args.email
                ]):
                    print(
                        f"Usage:\n{FDSA_BUILD_COMMAND_USAGE}\n\n"
                        "If no run config file is provided, the --species and --genome options must be set.\n"
                        "Use \"fdsa build --help\" for more information on usage."
                    )
                    exit()

            if build_args.run_config_path is not None:
                # Load run config if provided
                if not os.path.exists(build_args.run_config_path):
                    print(INVALID_CONFIG_PATH_MESSAGE)
                    exit()
                try:
                    run_config = ProgramRunConfig(
                        build_args.run_config_path,
                        internal_config,
                        user_config
                    )
                except ValueError as e:
                    print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
                    exit()

                # Allow args to override run config if provided
                threads = build_args.threads if build_args.threads else run_config.n_threads
                species = build_args.species if build_args.species else run_config.species
                genome = build_args.genome if build_args.genome else run_config.genome

            else:
                # Have previously validated that the first three are not None if build_args.run_config_path is not set
                species = build_args.species
                genome = build_args.genome
                threads = build_args.threads if build_args.threads else 1  # Default to single-threaded

            if build_args.email is None and user_config is None:
                print(NO_EMAIL_ADDRESS_MESSAGE)
                exit()

            email = build_args.email if build_args.email is not None else user_config.email_for_apis
            # Possible TODO: Final check for valid email_for_apis?

            if species not in internal_config.allowed_species:
                print(f"Invalid species (\"{species}\") specified. Options are: "
                      f"{'/'.join(internal_config.allowed_species)}")
                exit()

            if not os.path.exists(genome):
                print(f"The specified reference genome GTF file was not found. Path: {genome}")

            species_specific_data_dir = os.path.join(base_dir, "data", str(species))

            # Execute build
            _build(
                species_specific_data_dir,
                species,
                genome,
                email,
                threads,
                internal_config,
                user_config,
                force_redo_annotate_transcript_library=build_args.redo_annotation,
                force_regenerate_whole_genome_lookup=build_args.regenerate_lookup,
                force_regenerate_transcript_library=build_args.regenerate_transcripts,
                force_regenerate_genbank_features=build_args.regenerate_features,
                restart_partial_download=build_args.regenerate_partial_features_download
            )

        elif mode_arg.mode in ("run", "report", "filter"):

            # TODO: Implement --no-report flag

            enable_generate_report = False
            enable_generate_filtered_bam_files = False

            if mode_arg.mode == "run":
                parser = get_run_parser()
                usage = FDSA_RUN_COMMAND_USAGE
            elif mode_arg.mode == "report":
                enable_generate_report = True
                parser = get_report_parser()
                usage = FDSA_REPORT_COMMAND_USAGE
            else:
                enable_generate_filtered_bam_files = True
                parser = get_filter_parser()
                usage = FDSA_FILTER_COMMAND_USAGE

            if any([arg in subsequent_args for arg in ("-h", "--help")]):
                parser.print_help()
                exit()

            args = parser.parse_args(subsequent_args)

            if args.run_config_path is None:
                print(
                    f"usage: {usage}\n\n"
                    "Use \"fdsa run --help\" for more information on usage."
                )
                exit()

            run_config_path = args.run_config_path

            # Load run config
            try:
                run_config = ProgramRunConfig(
                    run_config_path,
                    internal_config,
                    user_config
                )
                if mode_arg.mode == "run":
                    if args.report:
                        enable_generate_report = True
                        run_config.generate_report = True
                        run_config.check_feature_counts()  # Ensure either featureCountsExecutable or geneCounts is set
                    if args.filter:
                        enable_generate_filtered_bam_files = True
                        run_config.generate_filtered_bam_files = True
                        run_config.check_samtools()  # Ensure samtoolsExecutable is set
            except ValueError as e:
                print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
                exit()

            if mode_arg.mode == "run":

                # Enable report generation if specified in config file and not via command line flag
                if run_config.generate_report:
                    enable_generate_report = True

                # The --no-report flag will override behaviour specified elsewhere
                if args.no_report:
                    enable_generate_report = False

                # Enable BAM file filtering if specified in config file and not via command line flag
                if run_config.generate_filtered_bam_files:
                    enable_generate_filtered_bam_files = True

                # The --no-filter flag will override behaviour specified elsewhere
                if args.no_filter:
                    enable_generate_filtered_bam_files = False

                species_specific_data_dir = os.path.join(base_dir, "data", run_config.species)

                # Execute run
                _run(
                    species_specific_data_dir,
                    internal_config,
                    run_config
                )

            if enable_generate_report:

                # Execute report
                _report(run_config)

            if enable_generate_filtered_bam_files:

                # Execute filter
                _filter(run_config)

        elif mode_arg.mode == "quant":

            quant_parser = get_quant_parser()
            quant_args = quant_parser.parse_args(subsequent_args)

            if quant_args.run_config_path is None and any([
                quant_args.genome is None,
                quant_args.input is None,
                quant_args.output is None
            ]):
                print(
                    f"Usage:\n{FDSA_QUANT_COMMAND_USAGE}\n\n"
                    "Use \"fdsa quant --help\" for more information on usage."
                )
                exit()

            # Load run config if path was supplied
            if quant_args.run_config_path is not None:
                try:
                    run_config = ProgramRunConfig(
                        quant_args.run_config_path,
                        internal_config,
                        user_config
                    )
                except ValueError as e:
                    print(CONFIG_FILE_FORMATTING_ERROR_MESSAGE + "\n" + str(e))
                    exit()
            else:
                run_config = None

            # Execute quant
            _quant(
                run_config,
                user_config,
                quant_args.stringtie,
                quant_args.prep_de,
                quant_args.genome,
                quant_args.input,
                quant_args.output,
                quant_args.read_length,
                quant_args.fdsa_results,
                quant_args.disable_assign_reference_gene,
                quant_args.threads
            )


if __name__ == "__main__":

    main()
