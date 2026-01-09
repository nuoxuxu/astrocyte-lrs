#!/bin/env python
import gffutils
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="store", dest="input")
    parser.add_argument("--output", action="store", dest="output")
    params = parser.parse_args()

    db = gffutils.create.create_db(
        data=params.input,
        dbfn=params.output,
        force=True,
        verbose=False,
        force_gff=False,
        checklines=1000,
        disable_infer_transcripts=True,
    )

if __name__ == "__main__":
    main()