
from argparse import ArgumentParser


def get_mode_parser(program_description: str) -> ArgumentParser:

    mode_parser = ArgumentParser(
        description=program_description,
        add_help=False
    )

    mode_parser.add_argument(
        "mode",
        type=str,
        choices=["user", "build", "run", "report"],
        help="\n".join([
            "user - (Required only once) Sets the user config file containing build information and optional default "
            "parameters.",
            "build - (Required once per species, before first run) Parses transcripts from the provided reference "
            "genome GTF, downloads feature annotation data, and links the two to enable directed splice analysis.",
            "run - Performs splice analysis according to parameters set in the specified run config file.",
            "report - Generates a graphical report from the output of \"fase run\"."
        ])
    )

    return mode_parser


def get_build_parser() -> ArgumentParser:

    build_parser = ArgumentParser()

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
        help="Species for which transcript annotation data will be downloaded. Must match the species of the "
             "provided reference genome GTF."
    )
    build_parser.add_argument(
        "-g",
        "--genome",
        type=str,
        required=False,
        help="Reference genome GTF from which transcripts will be parsed."
    )
    build_parser.add_argument(
        "-e",
        "--email",
        type=str,
        required=False,
        help="An identifying email address to be attached to GenBank API queries. Note: no emails will be "
             "sent - used for identification/anti-spam purposes only."
    )
    build_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="Number of threads to use during build process for speed improvements"
    )
    build_parser.add_argument(
        "--regenerate-lookup",
        action="store_true",
        help="Re-downloads transcript information from BioMart and regenerates the locally stored "
             "lookup used to convert between transcript IDs, gene IDs and gene names."
    )
    build_parser.add_argument(
        "--regenerate-transcripts",
        action="store_true",
        help="Re-parses transcripts from the provided reference genome GTF to create a stored transcript "
             "library."
    )
    build_parser.add_argument(
        "--regenerate-features",
        action="store_true",
        help="Re-downloads transcript feature annotation data from GenBank and stores the information locally."
    )
    build_parser.add_argument(
        "--redo-annotation",
        action="store_true",
        help="Re-applies feature annotation to the stored transcript library."
    )
    build_parser.add_argument(
        "--regenerate-all",
        action="store_true",
        help="Shorthand for --regenerate-lookup --regenerate-transcripts "
             "--regenerate-features --redo-annotation "
    )
    build_parser.add_argument(
        "--regenerate-partial-features-download",
        action="store_true",
        help="When combined with the --regenerate-features flag, forces any data partially downloaded data from "
             "GenBank to be discarded, restarting the download from scratch."
    )

    return build_parser


def get_run_parser() -> ArgumentParser:

    run_parser = ArgumentParser()

    run_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled below
        type=str,
        help=(
            "Path to a run config file containing a [RUN] section with the mandatory \"name\" (run name), \"feature\", "
            "\"species\", \"input\", \"bamEnding\", and \"output\" parameters set."
        )
    )
    run_parser.add_argument(
        "-r",
        "--report",
        action="store_true",
        help="Generates a visual report for the output of FASE's analysis upon completion."
    )

    return run_parser


def get_report_parser() -> ArgumentParser:

    run_parser = ArgumentParser()

    run_parser.add_argument(
        "run_config_path",
        nargs="?",
        default=None,  # Allow no value so missing value error can be manually handled below
        type=str,
        help=(
            "Path to a run config file containing a valid [RUN] section (see \"fase run --help\") and optional "
            "[REPORT], [SAMPLES], and [COLORS] sections."
            ""
        )
    )

    return run_parser
