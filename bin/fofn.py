#!/bin/env python
from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='Create a FOFN (file of file names) for each unique primer pair based on the file name.'
    )
    parser.add_argument('file_suffix', help='Suffix of the files to look for (e.g., ".flnc.filter_summary.report.json")')
    args = parser.parse_args()
    flist = list(Path("./").glob(args.file_suffix))

    def get_primer_pair(file_path):
        parts = file_path.name.split(".")
        primer_part = parts[1]
        return primer_part
    unique_primer_Pair = list(set([get_primer_pair(file_path) for file_path in flist]))
    my_dict = {
            primer: [str(file_path) for file_path in flist if get_primer_pair(file_path) == primer] 
            for primer in unique_primer_Pair
        }

    for primer in my_dict.keys():
        with open(f"{primer}.fofn", "w") as f:
            for file_path in my_dict[primer]:
                f.write(f"{file_path}\n")

if __name__ == "__main__":
    main()