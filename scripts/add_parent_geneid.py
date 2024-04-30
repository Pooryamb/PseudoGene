#!/usr/bin/python3
from gff_cols import gff_cols

def add_parent_geneid(gff):
    """This function assumes that parent lines of gff contain one of the following types:
    "protein_coding_gene","pseudogene", "ncRNA_gene"
    It also assumes that each parent line is followed by its related lines. Finally, it 
    extracts the id of the parent line in gff"""
    gff["isparent"] = gff["type"].isin(["protein_coding_gene","pseudogene", "ncRNA_gene"]).astype(int)
    gff["parent_id"] = gff["isparent"].cumsum()
    parent_col_values = ["protein_coding_gene","pseudogene", "ncRNA_gene"]
    parent_rows = gff[gff["type"].isin(parent_col_values)]

    parent_rows["gene_id"] = parent_rows["attributes"].str.split(";",expand=True)[0].str.replace("ID=", "")
    parent_rows = parent_rows[["parent_id", "gene_id"]].reset_index(drop=True)
    gff = gff.merge(parent_rows, on ="parent_id")
    return gff[gff_cols + ["gene_id"]]
