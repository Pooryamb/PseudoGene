#!/usr/bin/python3

# MUSCLE did not support stop codons inside the protein sequence. To use MUSCLE, we replaced the stop codons with X
# After MUSCLE is done, we have to replace the Xs with stop codons again. This script does the task.

import glob
import os
from fasta2dict import fasta2dict

aligned_fasta_paths = glob.glob("../data/manual_checking/*/*/seq4dnds/prot_orths.afa")

for aligned_fasta_path in aligned_fasta_paths:
    pg_dir_path = os.path.dirname(aligned_fasta_path)
    x_locs_str = open(f"{pg_dir_path}/x_locs.txt").read().split("\t")
    x_locs = [int(x) for x in x_locs_str if x != ""]

    original_order = fasta2dict(aligned_fasta_path.replace(".afa", ".fasta")).keys()
    aligned_dict = fasta2dict(aligned_fasta_path)

    processed_aln_path = aligned_fasta_path.replace(".afa", "_4p2n.afa")
    with open(processed_aln_path, 'w') as processed_fle:
        for i, seq_id in enumerate(original_order):
            if i==0:
                aa_pos = 0
                corrected_aas = []
                for aa in aligned_dict[seq_id]:
                    if aa.isalnum():
                        aa_pos += 1
                        if aa== "X" and (aa_pos -1  not in x_locs):
                            corrected_aas.append("*")
                            continue
                    corrected_aas.append(aa)
                pg_seq = "".join(corrected_aas)
                processed_fle.write(f">{seq_id}\n{pg_seq}\n")
            else:
                processed_fle.write(f">{seq_id}\n{aligned_dict[seq_id]}\n")

