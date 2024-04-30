import pandas as pd
import sys
from select_first_name_in_header import select_first_name_in_header
from fasta2dict import fasta2dict

def find_genomicregion_start(cds_start, margin=900): 
    """This is for selecting the 900 nucleotides before the start of the genomic fragment"""
    minstart = 1
    genomicstart = max(cds_start - margin, cds_start - ((cds_start - minstart)//3)*3)
    return genomicstart

def find_genomicregion_end(cds_end, contiglen, margin=900): 
    """This is for selecting the 900 nucleotides after the stop of the genomic fragment"""
    genomicend = min(cds_end + margin, cds_end + ((contiglen - cds_end)//3)*3)
    return genomicend
 
def find_locations_cds_and_flanking(contiglen, simplified_gff):
    prot_loc = simplified_gff[(simplified_gff["type"]=="protein_coding_gene")]
    margin = 900
    proteins_contiglen = prot_loc.merge(contiglen, on = "contig")
    proteins_contiglen["frag_start"] = proteins_contiglen.apply(lambda row: find_genomicregion_start(row["start"]), axis=1)
    proteins_contiglen["frag_end"] = proteins_contiglen.apply(lambda row: find_genomicregion_end(row["end"], row["contig_len"]),axis=1)
    fragment_regions = proteins_contiglen.drop(columns= ["contig_len", "type"])[["gene_id", "contig", "start", "end", "frag_start", "frag_end", "strand"]]
    fragment_regions = fragment_regions.rename(columns={"start": "cds_start", "end":"cds_end"})
    return fragment_regions