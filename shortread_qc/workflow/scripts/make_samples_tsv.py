#!/usr/bin/env python3
"""
make_samples_tsv.py — build config/samples.tsv from directories of raw FASTQs.

Finds R1 files, pairs each with its R2, strips Illumina suffixes
(_S12_L001_R1_001.fastq.gz, _R1.fq.gz, .R1.fastq.gz, _1.fastq.gz ...) to get the
sample name, and groups multiple lanes of one sample into comma-separated lists.
Review the output before running: sample names become file names everywhere.

If the same sample name turns up in more than one folder (e.g. a re-sequenced
library in a later delivery), the files are merged into one sample by default
(more depth) and a WARNING lists them. Use --split-folders to keep each folder's
reads as a separate sample, named <sample>__<folder date or name>, so each
delivery can be checked on its own first.

Usage: make_samples_tsv.py [--split-folders] RAW_DIR [RAW_DIR ...] > config/samples.tsv
"""
import os
import re
import sys
from collections import defaultdict

R1_PAT = re.compile(r"^(?P<stem>.+?)(?P<sl>_S\d+)?(?P<lane>_L\d{3})?[._](?:R)?1(?P<tail>(_\d{3})?\.f(?:ast)?q(?:\.gz)?)$")


DATED = re.compile(r"^(20\d{2})[._-](\d{2})[._-](\d{2})")


def folder_tag(path, roots):
    """Short label for the top-level folder a file came from (its date if dated)."""
    ap = os.path.abspath(path)
    for r in roots:
        r = os.path.abspath(r)
        if ap.startswith(r + os.sep):
            rel = os.path.relpath(ap, r).split(os.sep)
            for part in [os.path.basename(r)] + rel[:-1]:
                m = DATED.match(part)
                if m:
                    return f"{m.group(1)}{m.group(2)}{m.group(3)}"
            return re.sub(r"[^A-Za-z0-9._-]", "_", rel[0] if len(rel) > 1 else os.path.basename(r))
    return "x"


def main():
    args = sys.argv[1:]
    split = "--split-folders" in args
    dirs = [a for a in args if a != "--split-folders"]
    if not dirs:
        sys.exit(__doc__)
    groups = defaultdict(list)
    for d in dirs:
        for root, _, files in os.walk(d):
            for f in sorted(files):
                m = R1_PAT.match(f)
                if not m:
                    continue
                r1 = os.path.abspath(os.path.join(root, f))
                # R2 = same name with the read digit swapped at the matched position
                r1_idx = m.end("lane") if m.group("lane") else (m.end("sl") if m.group("sl") else m.end("stem"))
                r2_name = f[:r1_idx] + f[r1_idx:].replace("1", "2", 1)
                r2 = os.path.join(root, r2_name)
                if not os.path.exists(r2):
                    print(f"# WARNING: no R2 for {r1} (looked for {r2_name})", file=sys.stderr)
                    continue
                name = re.sub(r"[^A-Za-z0-9._-]", "_", m.group("stem"))
                if split:
                    name = f"{name}__{folder_tag(r1, dirs)}"
                groups[name].append((m.group("lane") or "", r1, os.path.abspath(r2)))
    for name, files in sorted(groups.items()):
        folders = sorted({os.path.dirname(f[1]) for f in files})
        if len(folders) > 1:
            print(f"# WARNING: {name} has files in {len(folders)} folders and will be MERGED "
                  f"(re-run with --split-folders to keep them apart):", file=sys.stderr)
            for fo in folders:
                print(f"#   {fo}", file=sys.stderr)
    print("sample\tr1\tr2")
    for name in sorted(groups):
        lanes = sorted(groups[name])
        print(f"{name}\t{','.join(l[1] for l in lanes)}\t{','.join(l[2] for l in lanes)}")


if __name__ == "__main__":
    main()
