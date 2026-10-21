#!/usr/bin/env python3
"""Collect one sample's QC numbers (fastp, host removal, seqkit, adapter checks)
into a single-row TSV. aggregate_tsv.py stacks the rows into qc_summary.tsv."""
import argparse
import csv
import json


def read_kv(path):
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return {r[0]: r[1] for r in rows[1:] if len(r) >= 2}


def pct(a, b):
    return round(100.0 * a / b, 3) if b else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastp-json", required=True)
    ap.add_argument("--host-counts", required=True)
    ap.add_argument("--seqkit", required=True)
    ap.add_argument("--adapters-raw", required=True)
    ap.add_argument("--adapters-clean", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fp = json.load(open(a.fastp_json))
    summ = fp.get("summary", {})
    before = summ.get("before_filtering", {})
    after = summ.get("after_filtering", {})
    filt = fp.get("filtering_result", {})
    cut = fp.get("adapter_cutting", {})

    raw_reads = before.get("total_reads", 0)
    fastp_reads = after.get("total_reads", 0)

    host = read_kv(a.host_counts)
    host_pairs = int(host["input_pairs"])
    nonhuman_pairs = int(host["nonhuman_pairs"])
    human_pairs = host_pairs - nonhuman_pairs

    # seqkit stats -T -a over the clean R1/R2
    with open(a.seqkit) as fh:
        clean = list(csv.DictReader(fh, delimiter="\t"))
    clean_r1_reads = int(clean[0]["num_seqs"]) if clean else -1
    clean_r2_reads = int(clean[1]["num_seqs"]) if len(clean) > 1 else -1
    final_bases = sum(int(r["sum_len"]) for r in clean)
    final_q30 = round(sum(float(r["Q30(%)"]) for r in clean) / len(clean), 2) if clean else "NA"
    final_avg_len = round(sum(float(r["avg_len"]) for r in clean) / len(clean), 1) if clean else "NA"

    ar = read_kv(a.adapters_raw)
    ac = read_kv(a.adapters_clean)

    row = {
        "sample": a.sample,
        "instrument_guess": ar.get("instrument_guess", "NA"),
        "raw_pairs": raw_reads // 2,
        "raw_bases": before.get("total_bases", "NA"),
        "raw_q30_rate": before.get("q30_rate", "NA"),
        "raw_looks_already_trimmed": ar.get("looks_already_trimmed", "NA"),
        "adapter_family_raw": ar.get("dominant_adapter_family", "NA"),
        "adapter_r1_consensus_raw": ar.get("adapter_r1_consensus", "NA"),
        "adapter_r2_consensus_raw": ar.get("adapter_r2_consensus", "NA"),
        "fastp_r1_adapter": cut.get("read1_adapter_sequence", "NA"),
        "fastp_r2_adapter": cut.get("read2_adapter_sequence", "NA"),
        "pct_pairs_readthrough_raw": ar.get("pct_pairs_readthrough", "NA"),
        "pct_pairs_readthrough_clean": ac.get("pct_pairs_readthrough", "NA"),
        "pct_reads_adapter_signature_raw": ar.get("pct_reads_with_dominant_signature", "NA"),
        "pct_reads_adapter_signature_clean": ac.get("pct_reads_with_dominant_signature", "NA"),
        "pct_polyG_r1_raw": ar.get("pct_r1_polyG", "NA"),
        "pct_polyG_r1_clean": ac.get("pct_r1_polyG", "NA"),
        "fastp_pct_reads_adapter_trimmed": pct(cut.get("adapter_trimmed_reads", 0), raw_reads),
        "fastp_duplication_rate": fp.get("duplication", {}).get("rate", "NA"),
        "fastp_insert_peak": fp.get("insert_size", {}).get("peak", "NA"),
        "fastp_low_quality_reads": filt.get("low_quality_reads", "NA"),
        "fastp_too_many_N_reads": filt.get("too_many_N_reads", "NA"),
        "fastp_too_short_reads": filt.get("too_short_reads", "NA"),
        "fastp_pass_pairs": fastp_reads // 2,
        "pct_pairs_removed_fastp": pct(raw_reads - fastp_reads, raw_reads),
        "host_input_pairs": host_pairs,
        "human_pairs_removed": human_pairs,
        "pct_pairs_human": pct(human_pairs, host_pairs),
        "final_pairs": nonhuman_pairs,
        "final_bases": final_bases,
        "final_mean_read_len": final_avg_len,
        "final_q30_pct": final_q30,
        "pct_pairs_retained_overall": pct(nonhuman_pairs, raw_reads // 2),
        # Sanity check that the filter fired as counted: reads written == pairs counted in the BAM
        "counts_consistent": clean_r1_reads == clean_r2_reads == nonhuman_pairs,
        "fastp_version": summ.get("fastp_version", "NA"),
    }
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    main()
