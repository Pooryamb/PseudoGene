#!/usr/bin/python3
#NOTE: This script works only on positive frames.

import pandas as pd
from break_dmndali2frames import break_dmndali2frames
from break_ali2genomictypes import break_ali2genomictypes


start_stop_path = sys.argv[1]
start_stop = pd.read_csv(start_stop_path, sep="\t")

alignmentspath = sys.argv[2]
brokenalipath  = sys.argv[3]

alidf = pd.read_csv(alignmentspath, sep="\t")

frame_broken = []
for index,row in alidf.iterrows():
    if index%10000==0:
        print(f"Index up to {index//1000}k processed")
    frame_broken.extend(break_dmndali2frames(row))

frame_broken_df = pd.DataFrame(frame_broken, columns =["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "qlen", "qseq", "sseq"] )
frame_broken_df = frame_broken_df.merge(start_stop, left_on = "qseqid", right_on = "gene_id").drop(columns= ["gene_id"])

type_broken = []
for index,row in frame_broken_df.iterrows():
    if index%10000==0:
        print(f"Index up to {index//1000}k processed")
    type_broken.extend(break_ali2genomictypes(row))

type_broken_df = pd.DataFrame(type_broken, columns =["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "qlen", "qseq", "sseq", "region_type"])

type_broken_df.to_csv(brokenalipath, sep="\t", index=None)