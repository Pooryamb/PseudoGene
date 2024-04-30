#!/usr/bin/python3

import sys
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
from read_gff import read_gff
from simplify_gff import simplify_gff
from find_genes_to_exclude import find_overlapping_and_intronic_proteins
from select_first_name_in_header import select_first_name_in_header
from fasta2dict import fasta2dict
from find_locations_cds_and_flanking import find_locations_cds_and_flanking
from select_subseq import select_subseq
from dict2fasta import dict2fasta
from get_start_stop_codon_poses import get_poses_reversed_frame, get_poses_original_frame

gff_path = sys.argv[1]
genome_fasta_path = sys.argv[2]

gff = read_gff(gff_path)
simplified_gff = simplify_gff(read_gff(gff_path))

to_exclude = find_overlapping_and_intronic_proteins(gff, simplified_gff)

genome_dict = fasta2dict(genome_fasta_path, select_first_name_in_header)

contig2len_dict = {x:len(y) for (x,y) in genome_dict.items()}
contig2len_df = pd.DataFrame.from_dict(data=contig2len_dict, orient="index").reset_index().rename(columns={"index": "contig", 0: "contig_len"})

fragment_coords = find_locations_cds_and_flanking(contig2len_df, simplified_gff)
fragment_coords = fragment_coords[~fragment_coords["gene_id"].isin(to_exclude["gene_id"])]
fragment_coords.to_csv("../data/fragment_locations.tsv", sep="\t", index=None)


start_stop_ori_frame = get_poses_original_frame(fragment_coords)
start_stop_ori_frame.to_csv("../data/start_stop_pos_ori.tsv", sep="\t", index=None)
start_stop_rev_frame = get_poses_reversed_frame(fragment_coords)
start_stop_rev_frame.to_csv("../data/start_stop_pos_rev.tsv", sep="\t", index=None)


selected_regions = select_subseq(genome_dict, fragment_coords)
dict2fasta(selected_regions, "../data/cds_and_flanking.fasta")