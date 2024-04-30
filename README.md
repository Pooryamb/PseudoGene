
# Introduction

This repository contains scripts designed to identify genes where frame correction, based on their alignment with a homologous protein, leads to a significantly extended protein length annotated with domains. This analysis can help with identification of pseudogenes mislabeled as functional genes.

## Prerequisites

We have successfully run the script using the following tools (version has been specified too):

- Python==3.11.7
- Diamond==2.1.8
- InterProScan==5.63-95.0
- Pandas==2.1.4
- HyPhy==2.5.2
- PhyML==20120412
- PAL2NAL==v14
- MUSCLE=v5
- GNU Parallel

Ensure that the path to these tools is specified in the `scripts/run.sh` file.

## Setup

Before executing the scripts, you must provide several input files:

- A GFF file containing coordinates of genomic fragments.
- A genome FASTA file.
- A FASTA file of annotated proteins.
- A FASTA file of proteins from closely related species (used for selection analysis).
- A FASTA file of genomic sequences from closely related species (used for selection analysis).
- A table summarizing the orthologs for each gene in the query organism relative to its close species.

Place all these files in the `data` directory and ensure their paths are correctly specified in the `scripts/run.sh` file. The current repository comes with the same data for _Trypanosoma_ _brucei_.

## Output

After processing, the scripts will generate directories within `./data/manual_checking` for genes that meet the filtering criteria. These directories are categorized by the alignment frame of the genomic fragments relative to the subject:

- Positive frame (`pf`): Query aligns with the subject with an ORF whose direction is similar to the previously predicted ORF of query.
- Negative frame (`nf`): Query aligns with the subject with an ORF whose direction is opposite to the previously predicted ORF of query.

Each gene-specific directory will contain the following files:

- `ori_f_doms.tsv`: Domains predicted on the original ORF of the query.
- `nonbroken_ali.tsv`: Alignment of the query with the subject (Diamond output).
- `broken_ali.tsv`: Alignment of the query with the subject showing each frame separately. If the corrected ORF includes flanking regions, each region (3' or 5') is displayed on a separate line.
- `subject_doms.tsv`: Domains predicted on the subject used for correcting the ORF of the query.
- `correctedseq_doms.tsv`: Domains predicted on the corrected ORF of the query.
- `orths.nal.RELAX.json`: Output from the selection analysis.

These files facilitate a detailed examination of the gene corrections and are available for further analysis by the user.
