
from typing import Dict
import os
import subprocess


def quantify_isoforms(
    stringtie_executable_path: str,
    prep_de_script_path: str,
    gtf_path: str,
    filtered_bam_dir: str,
    output_dir: str,
    bam_file_ending: str = ".bam",
    threads: int = 1,
    check_exit_code: bool = False
) -> None:

    if not os.path.isfile(gtf_path):
        raise ValueError(f"Invalid path supplied for reference genome: {gtf_path}")

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

    individual_stringtie_output_locations: Dict[str, str] = {}
    _n = 0
    for bam_path in bam_abs_paths:

        print(os.path.basename(bam_path))

        # output_filename = os.path.basename(bam_path).replace(bam_file_ending, ".gtf")
        _basename = os.path.basename(bam_path)
        output_filename = _basename[:_basename.index(".")] + ".gtf"

        individual_stringtie_output_locations[bam_path] = os.path.join(output_base_dir, "individual", output_filename)

        # stringtie -G GENOME -o INDIVIDUAL.gtf -p THREADS BAM
        subprocess.run(
            [
                stringtie_executable_path,
                "-G", gtf_path,
                "-o", individual_stringtie_output_locations[bam_path],
                "-p", str(threads),
                bam_path
            ],
            cwd=output_base_dir,
            encoding="utf8",
            check=check_exit_code
        )

        # DEBUG: Break after second iteration
        _n += 1
        if _n >= 2:
            break

    print("[stringtie] ... done")

    # Run stringtie in --merge mode to get the combined transcript GTF

    print("[stringtie] Merging results...")

    # Write merge list for --merge mode
    merge_list_path = os.path.join(output_base_dir, "merge_list.txt")
    with open(merge_list_path, "w") as f:
        f.write("\n".join(list(individual_stringtie_output_locations.values())))

    # stringtie --merge -G GENOME -o MERGED.gtf MERGE_LIST
    subprocess.run(
        [
            stringtie_executable_path,
            "--merge",
            "-G", gtf_path,
            "-o", os.path.join(output_base_dir, "merged.gtf"),
            merge_list_path
        ],
        cwd=output_base_dir,
        encoding="utf8",
        check=check_exit_code
    )

    print("[stringtie] ... done")

    # Run stringtie in -e -b mode
    # In this step, the combined GTF from stringtie --merge is used with -G rather than the reference genome

    pass
