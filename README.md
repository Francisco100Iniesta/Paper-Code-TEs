# [PROJECT / PAPER SHORT TITLE]

Scripts and analysis workflows associated with the manuscript:

**"[FULL MANUSCRIPT TITLE]"**

Authors: 
Francisco Iniesta-Martinez1,†, María Llamas-López1,†, Esther Navarro-Manzano1-3, Alba Rodríguez-Ródenas3 , José Padilla1, Marina Fuentes-Custodio1, Francisco Abad-Navarro2-3, Javier Cuenca-Guardiola2, Carlos Bravo-Pérez1, Bruno Ramos-Molina4, Lidia Sánchez-Alcoholado4 , Benedicte Stavik6-7, Gareth J Sullivan8 Per Morten Sandset6-7, María Eugenia Chollet6, María Eugenia de la Morena-Barrio1, María Luisa Lozano1, Jesualdo Tomás Fernández-Breis1-2,Guiomar Perez-de Nanclares5,  Javier Corral1, Belén de la Morena-Barrio1   
1 Department of Haematology, Hospital Universitario Morales Meseguer, Centro Regional de Hemodonación, University of Murcia, IMIB Pascual Parrilla, CIBERER-ISCIII, 30003, Murcia, Spain
2 Department of Informatics and Systems, University of Murcia, CEIR Campus Mare Nostrum, IMIB-Arrixaca, Faculty of Computer Science, Murcia, Spain.
3 LongSeq Applications SL., Murcia, Spain 
4 Obesity, Diabetes and Metabolism laboratory, Biomedical Research Institute of Murcia (IMIB), Murcia, Spain.
5 (Epi)genetics of Rare Diseases, Biobizkaia Health Research Institute, Cruces University Hospital, Cruces-Barakaldo, Bizkaia, Spain
6 Institute of Clinical Medicine, Department of Haematology, Oslo University Hospital, Oslo, Norway
7 Research Institute of Internal Medicine, Oslo University Hospital. Oslo, Norway.
8 Department of Molecular Medicine, Oslo University Hospital, Oslo, Norway
† These authors contributed equally to this work

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
