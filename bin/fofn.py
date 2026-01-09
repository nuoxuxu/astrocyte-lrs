#!/bin/env python
import polars as pl
from pathlib import Path

def main():
    flnc_bams = list(Path("./").glob("*.flnc.bam"))

    def get_primer_pair(file_path):
        parts = file_path.name.split(".")
        primer_part = parts[1]
        return primer_part
    unique_primer_Pair = list(set([get_primer_pair(file_path) for file_path in flnc_bams]))
    my_dict = {
            primer: [str(file_path) for file_path in flnc_bams if get_primer_pair(file_path) == primer] 
            for primer in unique_primer_Pair
        }

    for primer in my_dict.keys():
        with open(f"{primer}.fofn", "w") as f:
            for file_path in my_dict[primer]:
                f.write(f"{file_path}\n")

if __name__ == "__main__":
    main()