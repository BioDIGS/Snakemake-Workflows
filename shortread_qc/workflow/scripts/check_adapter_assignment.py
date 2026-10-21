#!/usr/bin/env python3
"""
check_adapter_assignment.py — before fastp runs, confirm a sample's reads carry the
adapter family it was assigned (config fastp.adapter_sets).

Compares identify_adapters.py's per-sample result (dominant adapter signature and
de novo consensus) with the assigned R1/R2 sequences.

Verdicts:
  PASS          observed family = assigned family
  LOW_EVIDENCE  too little adapter read-through to call (long-insert library);
                passes, flagged, since overlap-based trimming still runs
  NO_ASSIGNMENT sample uses auto-detection only; passes
  MISMATCH      a different family dominates -> exits 1 (fastp for this sample does
                not run) unless --warn-only

Writes a one-row TSV used for results/adapters/adapter_audit.tsv.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from identify_adapters import SIGNATURES  # noqa: E402


def family_of(seq):
    """Which known signature family an adapter sequence belongs to."""
    seq = (seq or "").upper()
    if not seq:
        return ""
    for name, sig in SIGNATURES.items():
        n = min(len(sig), len(seq), 13)
        if seq[:n] == sig[:n]:
            return name
    return "unknown/custom"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--report", required=True, help="identify_adapters.py --tsv output")
    ap.add_argument("--set-label", default="")
    ap.add_argument("--r1", default="")
    ap.add_argument("--r2", default="")
    ap.add_argument("--min-evidence-pct", type=float, default=0.05,
                    help="min %% of reads carrying the dominant signature to make a call")
    ap.add_argument("--warn-only", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.report) as fh:
        kv = {r["key"]: r["value"] for r in csv.DictReader(fh, delimiter="\t")}
    observed = kv.get("dominant_adapter_family", "none")
    try:
        evidence = float(kv.get("pct_reads_with_dominant_signature", 0))
    except ValueError:
        evidence = 0.0
    consensus = kv.get("adapter_r1_consensus", "NA")
    assigned = family_of(a.r1)

    if not a.r1 and not a.r2:
        verdict, note = "NO_ASSIGNMENT", "auto-detection only"
    elif observed == "none" or evidence < a.min_evidence_pct:
        verdict = "LOW_EVIDENCE"
        note = (f"only {evidence:.3f}% of reads carry an adapter signature; cannot confirm "
                f"(overlap-based trimming still applies)")
    elif assigned == "unknown/custom":
        # custom adapter: accept if it agrees with the de novo consensus
        ok = consensus != "NA" and (a.r1[:12] in consensus or consensus[:12] in a.r1)
        verdict = "PASS" if ok else "MISMATCH"
        note = f"custom adapter vs consensus {consensus}"
    elif observed == assigned:
        verdict, note = "PASS", ""
    else:
        verdict = "MISMATCH"
        note = (f"assigned {a.set_label} ({assigned}) but reads show {observed} in {evidence:.2f}% of reads; "
                f"consensus {consensus}")

    row = {"sample": a.sample, "verdict": verdict, "assigned_set": a.set_label or "-",
           "assigned_family": assigned or "-", "observed_family": observed,
           "pct_reads_observed_signature": evidence, "observed_r1_consensus": consensus,
           "instrument_guess": kv.get("instrument_guess", "NA"), "note": note}
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader()
        w.writerow(row)
    msg = f"{a.sample}: {verdict}" + (f" — {note}" if note else "")
    print(msg, file=sys.stderr)
    if verdict == "MISMATCH" and not a.warn_only:
        sys.exit("Adapter MISMATCH: fix fastp.adapter_sets for this sample (or set "
                 "fastp.adapter_check_warn_only: true). fastp was not run.")


if __name__ == "__main__":
    main()
