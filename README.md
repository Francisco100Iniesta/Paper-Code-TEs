# [PROJECT / PAPER SHORT TITLE]

Scripts and analysis workflows associated with the manuscript:

**"[FULL MANUSCRIPT TITLE]"**

Authors: [Authors]  
Journal: Genome Biology  
Status: Manuscript in preparation / submitted

## Overview

This repository contains the scripts and computational workflows used to
generate the analyses presented in the manuscript.

The analyses include:

- whole-genome sequence processing
- DNA methylation analysis
- variant analysis
- genomic and epigenomic summary statistics
- statistical analyses
- generation of manuscript figures and tables

This repository contains analysis code only. Individual-level human genomic
data are not distributed through GitHub.

## Data availability

The study includes whole-genome sequencing and DNA methylation data from
12 human individuals.

Individual-level genomic data, including sequencing/alignment and
methylation-related files, are being deposited in the
European Genome-phenome Archive (EGA) under controlled access.

EGA study accession:

`EGASXXXXXXXXXXX`

[Replace with the final accession once available.]

Access to controlled data is subject to approval by the corresponding
Data Access Committee (DAC).

No personally identifiable or individual-level genomic data are included
in this repository.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── environment.yml
├── config/
│   └── config.yaml
│
├── scripts/
│   ├── 01_preprocessing/
│   ├── 02_methylation/
│   ├── 03_variants/
│   ├── 04_statistics/
│   └── 05_figures/
│
├── workflow/
│   └── [Snakefile / main.nf]
│
├── metadata/
│   └── example_metadata.tsv
│
└── results/
    └── README.md
