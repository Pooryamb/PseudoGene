#!/usr/bin/python3
import pandas as pd
from dmnd_header_fmt6 import dmnd_header_fmt6

dmnd_ali_prefix= "../data/genome_vs_proteome"
pairwise_ali_path = dmnd_ali_prefix + "_pairwise.txt"
tabular_ali_path = dmnd_ali_prefix + "_fmt6.txt"
tabular_and_alis_path = dmnd_ali_prefix + "_fmt6_alis.txt"

pw_seq_alis = []

with open(pairwise_ali_path) as file:
    content = file.read()

text_parts = content.split("Query= ")[1:]
for part in text_parts:
    lines_of_gene = part.split("\n")
    qseqid  = lines_of_gene[0].strip()
    hit_parts = part.split(">")[1:]
    
    for hit in hit_parts:
        qseq = ""
        sseq = ""

        sseqid  = hit.split("\n")[0].strip()
        seq_lines = hit.split("\n")[7:]
        seq_lines = [x for x in seq_lines if len(x.strip()) > 0]
    
        for line in seq_lines:
            if line.startswith("Query"):
                qseq = qseq + line.split()[2]
            elif line.startswith("Sbjct"):
                sseq = sseq + line.split()[2]
        pw_seq_alis.append([qseqid, sseqid, qseq, sseq])

pw_seq_alis_df = pd.DataFrame(pw_seq_alis, columns = ["qseqid", "sseqid", "qseq", "sseq"])
fmt6_df = pd.read_csv(tabular_ali_path, sep="\t", header=None, names = dmnd_header_fmt6)

fmt6_alnseq = pd.concat([fmt6_df,pw_seq_alis_df[["qseq", "sseq"]]], axis=1)
fmt6_alnseq = fmt6_alnseq.drop_duplicates(["qseqid", "sseqid"])

fmt6_alnseq.to_csv(tabular_and_alis_path, sep="\t", index=None)