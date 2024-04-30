#!/usr/bin/python3

import pandas as pd
import sys
from read_gff import read_gff
from add_parent_geneid import add_parent_geneid

def get_start_end(df):
    if df.iloc[0]["type"]== "protein_coding_gene":
        cds = df[df["type"]=="CDS"]
        start = cds["start"].min()
        end = cds["end"].max()
        return start,end
    start = df.iloc[0]["start"]
    end = df.iloc[0]["end"]
    return start,end

def fmt_genomic_frag(df):
    gene_id = df.iloc[0]["gene_id"]
    contig = df.iloc[0]["seqid"]
    strand=df.iloc[0]["strand"]
    start, end = get_start_end(df)
    type_region = df.iloc[0]["type"]
    return pd.Series({"gene_id":gene_id, "contig":contig, "start":start, "end":end, "strand":strand, "type":type_region})

def simplify_gff(gff_df):
    gff_parent = add_parent_geneid(gff_df)
    simplified_table = gff_parent.groupby("gene_id").apply(fmt_genomic_frag)
    simplified_table.reset_index(drop=True, inplace=True)
    return simplified_table
    