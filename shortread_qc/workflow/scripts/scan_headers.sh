#!/usr/bin/env bash
# scan_headers.sh FASTQ[.gz] OUT.tsv
# One pass over every read header. Counts reads per
# (instrument, run, flowcell, lane, filter flag, index) so the sample sheet can
# report exactly which runs/flowcells/lanes a sample came from and in what share.
set -euo pipefail
in=$1; out=$2
if [[ $in == *.gz ]]; then
  if command -v pigz &>/dev/null; then dec=(pigz -dc -p 2); else dec=(gzip -dc); fi
else
  dec=(cat)
fi
{
  printf "instrument\trun\tflowcell\tlane\tfilter\tindex\treads\n"
  "${dec[@]}" "$in" | awk '
    NR % 4 == 1 {
      n = split($1, a, ":")
      inst = substr(a[1], 2)
      if (n >= 7) { run = a[2]; fc = a[3]; lane = a[4] }
      else        { run = "NA"; fc = "NA"; lane = "NA" }
      flt = "NA"; idx = "NA"
      if (NF > 1 && $2 ~ /^[12]:[YN]:/) { m = split($2, b, ":"); flt = b[2]; idx = b[m] }
      if (n < 7) inst = "NON_ILLUMINA_HEADER"
      c[inst "\t" run "\t" fc "\t" lane "\t" flt "\t" idx]++
    }
    END { for (k in c) print k "\t" c[k] }' | sort
} > "$out"
