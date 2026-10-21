#!/usr/bin/env python3
"""
validate_spike.py — score the human spike-in check.

Counts read pairs by origin (simulated human reads carry TAG in the flowcell
field of the read name) in the pre-QC mix and in the final clean R1, then
reports how many human pairs leaked through and how many soil pairs were lost.

Soil-pair loss includes fastp filtering (quality/length) as well as any false
human assignments; the per-sample fastp numbers separate the two.
"""
import argparse
import gzip


def count(path, tag):
    tagged = other = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            name = line[1:].split()[0]
            f = name.split(":")
            if len(f) > 2 and f[2] == tag:
                tagged += 1
            else:
                other += 1
    return tagged, other


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mix-r1", required=True)
    ap.add_argument("--clean-r1", required=True)
    ap.add_argument("--tag", default="HSSPIKE")
    ap.add_argument("--max-leak-pct", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    spike_in, soil_in = count(a.mix_r1, a.tag)
    spike_out, soil_out = count(a.clean_r1, a.tag)
    leak = 100.0 * spike_out / spike_in if spike_in else float("nan")
    soil_kept = 100.0 * soil_out / soil_in if soil_in else float("nan")
    verdict = "PASS" if spike_in and leak <= a.max_leak_pct else "FAIL"

    lines = [
        "Human spike-in validation",
        "=" * 60,
        f"Simulated human pairs in mix     : {spike_in:,}",
        f"Simulated human pairs in output  : {spike_out:,}  ({leak:.3f}% leaked)",
        f"Real soil pairs in mix           : {soil_in:,}",
        f"Real soil pairs in output        : {soil_out:,}  ({soil_kept:.2f}% retained;",
        "                                    loss = fastp filtering + any false human hits)",
        f"Threshold                        : leak <= {a.max_leak_pct}%",
        f"VERDICT                          : {verdict}",
    ]
    text = "\n".join(lines)
    print(text)
    with open(a.out, "w") as fh:
        fh.write(text + "\n")


if __name__ == "__main__":
    main()
