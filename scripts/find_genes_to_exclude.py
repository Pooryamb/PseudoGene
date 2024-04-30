#!/usr/bin/python3

import pandas as pd
import sys
import os
from read_gff import read_gff
from add_parent_geneid import add_parent_geneid
from overlap import overlap

def find_overlapping_proteins(locations):
    locations_x_y = locations.merge(locations, on = "contig")
    locations_x_y = locations_x_y[(locations_x_y["type_x"]== "protein_coding_gene") & (locations_x_y["gene_id_x"] != locations_x_y["gene_id_y"])]
    locations_x_y["overlap"] = overlap(locations_x_y, ["start_x", "end_x"], ["start_y", "end_y"])
    overlappings = locations_x_y[locations_x_y["overlap"]>0][["gene_id_x"]].rename(columns={"gene_id_x":"gene_id"}).drop_duplicates().reset_index(drop=True)
    return overlappings
    
def contains_intron(gff_cds):
    """This function takes the CDS lines of each gene_id as input and outputs if the gene has multiple exons"""
    cds_cds = gff_cds.merge(gff_cds, on = "type")
    cds_cds["overlap"] = overlap(cds_cds, ["start_x", "end_x"], ["start_y", "end_y"])
    return (cds_cds["overlap"]<=0).any()
    
def find_intron_containing_proteins(gff):
    """The input of this function should be a gff file + a column describing the gene id of each row. The gene_id can be added by 
    add_parent function"""
    gff_cds = gff[gff["type"]=="CDS"]
    genes_intron_info = gff_cds.groupby("gene_id").apply(contains_intron).reset_index()
    genes_intron_info.columns = ["gene_id", "has_intron"]
    return genes_intron_info[genes_intron_info["has_intron"]][["gene_id"]].reset_index(drop=True)
    
def find_overlapping_and_intronic_proteins(gff, simplified_gff):
    gff = add_parent_geneid(gff)
    overlapping_proteins = find_overlapping_proteins(simplified_gff)
    intron_containing_proteins = find_intron_containing_proteins(gff)
    to_exclude = pd.concat([overlapping_proteins, intron_containing_proteins]).reset_index(drop=True).drop_duplicates()
    return to_exclude
