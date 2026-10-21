# BioDIGS soil metagenomics workflows

The Snakemake workflows used to process and analyze sequencing data for the BioDIGS soil metagenomics project. Each workflow is in its own folder with its own README, configuration and pinned conda environments, and can be run on its own.

## Workflows

| Workflow | What it does |
|---|---|
| [`shortread_qc/`](shortread_qc/) | Illumina read cleaning: adapter and poly-G trimming, quality filtering, human-read removal against T2T-CHM13 |
| [`hifi_assembly/`](hifi_assembly/) | Per-sample PacBio HiFi assembly with metaMDBG and myloasm, with read and contig QC |
| [`sylph/`](sylph/) | Read-level taxonomic profiling of short and long reads with sylph |

## How they connect

```
Illumina reads ──► shortread_qc ──► clean reads ──┐
                                                   ├──► sylph
HiFi reads ──► human removal ──► clean reads ─────┤
                                                   └──► hifi_assembly ──► contigs
```

## Reproducing results

Each workflow's README gives the inputs, reference databases and commands needed to rerun it. Every run records the tool versions, parameters, configuration and input/output checksums it used in that workflow's provenance folder.

## Data availability

Raw and cleaned reads: *[accessions to be added]*

## Citation and license

*[To be added]*
