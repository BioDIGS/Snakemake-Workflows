#!/usr/bin/env python3
"""Stack TSVs that share a header: aggregate_tsv.py OUT IN1 [IN2 ...]"""
import sys

out, ins = sys.argv[1], sys.argv[2:]
header = None
with open(out, "w") as o:
    for p in ins:
        with open(p) as fh:
            lines = fh.read().splitlines()
        if not lines:
            continue
        if header is None:
            header = lines[0]
            o.write(header + "\n")
        elif lines[0] != header:
            sys.exit(f"ERROR: header mismatch in {p}")
        for ln in lines[1:]:
            o.write(ln + "\n")
