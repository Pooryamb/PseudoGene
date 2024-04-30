#!/usr/bin/python3


import sys
import os
import glob
import sys
from pathlib import Path
import pandas as pd
from fasta2dict import fasta2dict
from translate import translate
from revcomp import revcomp


candidates_path = glob.glob("../data/manual_checking/*/*")

candidates_id = [os.path.basename(x) for x in candidates_path]


genus_proteome_path = sys.argv[1]
genus_cds_path =  sys.argv[2]
orth_path = sys.argv[3]

prot_fasta_dict = fasta2dict(genus_proteome_path)
cds_fasta_dict = fasta2dict(genus_cds_path)
cds_and_flanking_fasta_dict = fasta2dict("../data/cds_and_flanking.fasta")
cds_and_flanking_revcom_fasta_dict = {x: revcomp(y) for x,y in cds_and_flanking_fasta_dict.items()}

def get_all_ids_of_orths(ortho_list):
    non_null_orthologs = ortho_list[~ortho_list.isnull()]
    ortho_list_single = []
    for gene in non_null_orthologs:
        ortho_list_single.extend(gene.split(", "))
    return ortho_list_single

orth_data_trypanosomes = pd.read_csv(orth_path, sep="\t")
for path in candidates_path:
    pg, functional_copy = open(f"{path}/nonbroken_ali.tsv").readlines()[1].split("\t")[:2]
    peers = orth_data_trypanosomes[orth_data_trypanosomes[orth_data_trypanosomes.columns[0]]==functional_copy]
    if peers.shape[0] == 0:
        continue
    peers_list = get_all_ids_of_orths(peers.iloc[0])
    seq_dir_path = f"{path}/seq4dnds"
    Path(seq_dir_path).mkdir(parents=True, exist_ok=True)
    if "/pf/" in path:
        cds_flanking = cds_and_flanking_fasta_dict
        broken_filename = "broken_ali.tsv"
    else:
        cds_flanking = cds_and_flanking_revcom_fasta_dict
        broken_filename = "broken_rev_ali.tsv"
        
    with open(f"{seq_dir_path}/cds_orths.fasta" ,'w') as cds_orths_fle, \
         open(f"{seq_dir_path}/prot_orths.fasta" ,'w') as prot_orths_fle:

        broken_ali = pd.read_csv(f"{path}/{broken_filename}", sep="\t")

        cds_seq = ""
        for index, row in broken_ali.iterrows():
            start, end = row["qstart"], row["qend"]
            cds_seq = cds_seq + cds_flanking[pg][start-1: end]

        cds_orths_fle.write(f">{pg}\n{cds_seq}\n")
        prot_seq = translate(cds_seq)
        x_indices = [str(i) for i, x in enumerate(prot_seq) if x == "X"]
        with open(f"{seq_dir_path}/x_locs.txt" ,'w') as x_locs_fle:
            x_locs_fle.write("\t".join(x_indices))
            
        prot_seq = prot_seq.replace("*", "X")
        prot_orths_fle.write(f">{pg}\n{prot_seq}\n")

        for gene_id in peers_list:
            cds_orths_fle.write(f">{gene_id}\n{cds_fasta_dict[gene_id]}\n")
            prot_orths_fle.write(f">{gene_id}\n{prot_fasta_dict[gene_id]}\n")