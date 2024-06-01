
from typing import Dict, Optional
import os
import subprocess

from src.analysis.core import name_output
from src.downstream.process_stringtie import format_stringtie_matrices, annotate_formatted_stringtie_results


def quantify_isoforms(
    stringtie_executable_path: str,
    prep_de_script_path: str,
    reference_gtf_path: str,
    filtered_bam_dir: str,
    output_base_dir: str,
    bam_file_ending: str = ".bam",
    read_length: Optional[int] = None,
    fase_results_path: Optional[str] = None,
    assign_reference_gene: bool = True,
    threads: int = 1,
    check_exit_code: bool = False
) -> None:

    # TODO: Add option to delete intermediary files & directories upon completion

    if not os.path.isfile(reference_gtf_path):
        raise ValueError(f"Invalid path supplied for reference genome: {reference_gtf_path}")

    if not os.path.isdir(filtered_bam_dir):
        raise ValueError(f"Invalid path supplied for filtered BAM directory: {filtered_bam_dir}")

    stringtie_out_dir = os.path.join(output_base_dir, "STRINGTIE")
    if not os.path.isdir(stringtie_out_dir):
        if os.path.exists(stringtie_out_dir):
            raise FileExistsError(f"Could not create output directory ({stringtie_out_dir}): already exists as file")
        os.makedirs(stringtie_out_dir)

    _bam_file_ending_length = len(bam_file_ending)
    bam_abs_paths = [
        os.path.join(os.path.abspath(filtered_bam_dir), f)
        for f in os.listdir(filtered_bam_dir)
        if f[-_bam_file_ending_length:] == bam_file_ending
    ]
    if len(bam_abs_paths) == 0:
        raise ValueError(f"No BAM files found in supplied input directory: {filtered_bam_dir}")

    # Run stringtie individually for each sample

    print("[Stringtie] Assembling transcripts from individual BAM files...")

    assembly_dir = os.path.join(stringtie_out_dir, "assembly")
    stringtie_assembly_output_locations: Dict[str, str] = {}

    _n = 0
    for bam_path in bam_abs_paths:

        print(os.path.basename(bam_path))

        _basename = os.path.basename(bam_path)
        _output_local_filename = _basename[:_basename.index(".")] + ".gtf"

        stringtie_assembly_output_locations[bam_path] = os.path.join(assembly_dir, _output_local_filename)

        # stringtie -G GENOME -o INDIVIDUAL.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-G", reference_gtf_path,
                "-o", stringtie_assembly_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=stringtie_out_dir,
            encoding="utf8",
            check=check_exit_code
        )

        # DEV: Break after second iteration
        _n += 1
        if _n >= 2:
            print("[DEV] stopped after 2 files")
            break

    print("[Stringtie] ... done")

    # Run stringtie in --merge mode to get the combined transcript GTF

    print("[Stringtie] Merging results...")

    # Write merge list for --merge mode
    merge_list_path = os.path.join(assembly_dir, "merge_list.txt")
    with open(merge_list_path, "w") as f:
        f.write("\n".join(list(stringtie_assembly_output_locations.values())))

    merged_gtf_path = os.path.join(stringtie_out_dir, "merged.gtf")

    # stringtie --merge -G GENOME -o MERGED.gtf MERGE_LIST
    subprocess.run(
        [
            stringtie_executable_path,
            "--merge",
            "-G", reference_gtf_path,
            "-o", merged_gtf_path,
            merge_list_path
        ],
        cwd=stringtie_out_dir,
        encoding="utf8",
        check=check_exit_code
    )

    print("[Stringtie] ... done")

    # Run stringtie in -e -b mode
    # In this step, the combined GTF from stringtie --merge is used with -G rather than the reference genome

    print("[Stringtie] Quantifying transcripts...")

    quantified_dir = os.path.join(stringtie_out_dir, "quantified")
    stringtie_quantified_output_locations: Dict[str, str] = {}

    _n = 0
    for bam_path in bam_abs_paths:

        _basename = os.path.basename(bam_path)
        _basename_stripped = _basename[:_basename.index(".")]
        _output_local_filename = _basename_stripped + ".quant.gtf"

        stringtie_quantified_output_locations[bam_path] = os.path.join(
            quantified_dir, _basename_stripped, _output_local_filename
        )

        print(_basename)

        # stringtie -e -B -G MERGED.gtf -o OUT_QUANT.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-e", "-B",
                "-G", merged_gtf_path,
                "-o", stringtie_quantified_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=stringtie_out_dir,
            encoding="utf8",
            check=check_exit_code
        )

        # DEV: Break after second iteration
        _n += 1
        if _n >= 2:
            print("[DEV] stopped after 2 files")
            break

    print("[Stringtie] ... done")

    # Run prepDE script

    print("[Stringtie] Generating transcript count matrix (prepDE.py3)...")

    prep_de_out_dir = os.path.join(stringtie_out_dir, "prepDE")

    if not os.path.isdir(prep_de_out_dir):
        if os.path.isfile(prep_de_out_dir):
            raise ValueError(f"Could not create prepDE output directory {prep_de_out_dir} - File already exists")
        os.makedirs(prep_de_out_dir)

    prep_de_gene_counts_path = os.path.join(prep_de_out_dir, "gene_count_matrix.csv")
    prep_de_transcript_counts_path = os.path.join(prep_de_out_dir, "transcript_count_matrix.csv")

    # prepDE.py3 -i INPUT [-g GENE_COUNTS.csv] [-t TRANSCRIPT_COUNTS.csv] [-l READ_LENGTH]
    prep_de_cmd_split = [
        "python3", prep_de_script_path,
        "-i", quantified_dir,
        "-g", prep_de_gene_counts_path,
        "-t", prep_de_transcript_counts_path
    ]
    if read_length is not None:
        prep_de_cmd_split += ["-l", read_length]
    subprocess.run(
        prep_de_cmd_split,
        cwd=prep_de_out_dir,
        encoding="utf8",
        check=check_exit_code
    )

    print("[Stringtie] ... done")

    print("[FASE] Combining stringtie results...")

    formatted_stringtie_output_path = os.path.join(stringtie_out_dir, "combined_stringtie_results.csv")

    format_stringtie_matrices(
        prep_de_gene_counts_path,
        prep_de_transcript_counts_path,
        merged_gtf_path,
        formatted_stringtie_output_path
    )

    print("[FASE] ... done")

    if fase_results_path is not None:

        print("[FASE] Annotating stringtie results with feature overlap...")

        annotate_formatted_stringtie_results(
            formatted_stringtie_output_path,
            fase_results_path,
            assign_reference_gene=assign_reference_gene
        )

        print("[FASE] ... done")

    else:

        print("[FASE] No FASE output supplied - stringtie results will not be annotated.")
