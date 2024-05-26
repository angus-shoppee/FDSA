
from typing import Dict
import os
import subprocess


def quantify_isoforms(
    stringtie_executable_path: str,
    prep_de_script_path: str,
    reference_gtf_path: str,
    filtered_bam_dir: str,
    output_dir: str,
    bam_file_ending: str = ".bam",
    threads: int = 1,
    check_exit_code: bool = False
) -> None:

    if not os.path.isfile(reference_gtf_path):
        raise ValueError(f"Invalid path supplied for reference genome: {reference_gtf_path}")

    if not os.path.isdir(filtered_bam_dir):
        raise ValueError(f"Invalid path supplied for filtered BAM directory: {filtered_bam_dir}")

    output_base_dir = os.path.join(output_dir, "STRINGTIE")
    if not os.path.isdir(output_base_dir):
        if os.path.exists(output_base_dir):
            raise FileExistsError(f"Could not create output directory ({output_base_dir}): already exists as file")
        os.makedirs(output_base_dir)

    _bam_file_ending_length = len(bam_file_ending)
    bam_abs_paths = [
        os.path.join(os.path.abspath(filtered_bam_dir), f)
        for f in os.listdir(filtered_bam_dir)
        if f[-_bam_file_ending_length:] == bam_file_ending
    ]
    if len(bam_abs_paths) == 0:
        raise ValueError(f"No BAM files found in supplied input directory: {filtered_bam_dir}")

    # Run stringtie individually for each sample

    print("[stringtie] Assembling transcripts from individual BAM files...")

    stringtie_assembly_output_locations: Dict[str, str] = {}
    _n = 0
    for bam_path in bam_abs_paths:

        print(os.path.basename(bam_path))

        _basename = os.path.basename(bam_path)
        _output_local_filename = _basename[:_basename.index(".")] + ".gtf"

        stringtie_assembly_output_locations[bam_path] = os.path.join(
            output_base_dir, "assembly", _output_local_filename
        )

        # stringtie -G GENOME -o INDIVIDUAL.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-G", reference_gtf_path,
                "-o", stringtie_assembly_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=output_base_dir,
            encoding="utf8",
            check=check_exit_code
        )

        # DEV: Break after second iteration
        _n += 1
        if _n >= 2:
            print("[DEV] stopped after 2 files")
            break

    print("[stringtie] ... done")

    # Run stringtie in --merge mode to get the combined transcript GTF

    print("[stringtie] Merging results...")

    # Write merge list for --merge mode
    merge_list_path = os.path.join(output_base_dir, "merge_list.txt")
    with open(merge_list_path, "w") as f:
        f.write("\n".join(list(stringtie_assembly_output_locations.values())))

    merged_gtf_path = os.path.join(output_base_dir, "merged.gtf")

    # stringtie --merge -G GENOME -o MERGED.gtf MERGE_LIST
    subprocess.run(
        [
            stringtie_executable_path,
            "--merge",
            "-G", reference_gtf_path,
            "-o", merged_gtf_path,
            merge_list_path
        ],
        cwd=output_base_dir,
        encoding="utf8",
        check=check_exit_code
    )

    print("[stringtie] ... done")

    # Run stringtie in -e -b mode
    # In this step, the combined GTF from stringtie --merge is used with -G rather than the reference genome

    print("[stringtie] Quantifying transcripts...")

    stringtie_assembly_output_locations: Dict[str, str] = {}

    _n = 0
    for bam_path in bam_abs_paths:

        _basename = os.path.basename(bam_path)
        _output_local_filename = _basename[:_basename.index(".")] + "_quant.gtf"

        stringtie_assembly_output_locations[bam_path] = os.path.join(
            output_base_dir, "quantified", _output_local_filename
        )

        print(_basename)

        # stringtie -e -B -G MERGED.gtf -o OUT_QUANT.gtf BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-e", "-B",
                "-G", merged_gtf_path,
                "-o", stringtie_assembly_output_locations[bam_path],
                bam_path
            ]
        )

        # DEV: Break after second iteration
        _n += 1
        if _n >= 2:
            print("[DEV] stopped after 2 files")
            break

    print("[stringtie] ... done")

    pass
