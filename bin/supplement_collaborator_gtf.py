#!/usr/bin/env python3

import argparse
import re


def get_transcript_id(line):
    m = re.search(r'transcript_id "([^"]+)"', line)
    return m.group(1) if m else None


def get_base_id(tid):
    return tid.rsplit('_', 1)[0] if '_' in tid else tid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("collaborator_gtf")
    parser.add_argument("final_gtf")
    parser.add_argument("-o", "--output", default="supplemented_collaborator.gtf")
    args = parser.parse_args()

    collab_lines = []
    collab_tids = set()
    with open(args.collaborator_gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            collab_lines.append(line)
            tid = get_transcript_id(line)
            if tid:
                collab_tids.add(tid)

    collab_base_ids = {get_base_id(t) for t in collab_tids}

    with open(args.output, "w") as out:
        for line in collab_lines:
            out.write(line)
        with open(args.final_gtf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tid = get_transcript_id(line)
                if tid is None or tid not in collab_base_ids:
                    out.write(line)

    # Assertions
    result_base_ids = set()
    result_underscore_tids = set()
    with open(args.output) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tid = get_transcript_id(line)
            if tid:
                result_base_ids.add(get_base_id(tid))
                if '_' in tid:
                    result_underscore_tids.add(tid)

    final_tids = set()
    with open(args.final_gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tid = get_transcript_id(line)
            if tid:
                final_tids.add(tid)

    assert result_base_ids == final_tids, (
        "Base ID set mismatch: "
        + str(len(result_base_ids - final_tids)) + " extra in result, "
        + str(len(final_tids - result_base_ids)) + " missing from result"
    )

    assert result_underscore_tids == collab_tids, (
        "Underscore ID set mismatch: "
        + str(len(result_underscore_tids - collab_tids)) + " extra in result, "
        + str(len(collab_tids - result_underscore_tids)) + " missing from result"
    )

    print("All assertions passed.")


if __name__ == "__main__":
    main()
