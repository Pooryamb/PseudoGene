#!/usr/bin/python3

import pandas as pd
def region_type(loc, start_codon, stop_codon):
    if loc < start_codon:
        return "five"
    if loc <=stop_codon:
        return "coding"
    return "three"
    
def break_ali2genomictypes(ali_start_stop): #Here, row must be the alignment row + start and stop codon positions
    type_broken = []
    for index, row in ali_start_stop.iterrows():
        start_codon, stop_codon = row["start_codon"], row["stop_codon"]
        partitioned_alis = []   #Each section will be represented by "qstart", "qend", "sstart", "send", "qframe", "qlen", "qseq", "sseq", "region_type"
        q_fr_start,s_fr_start,qframe = row["qstart"], row["sstart"], row["qframe"] #_fr_start shows the frame that the alignment starts with
        q_cur, s_cur = q_fr_start, s_fr_start
        last_region_type = region_type(q_cur, start_codon,stop_codon)

        #q_cur and s_cur is for keeping the record of the position of the cursor on query and subject
        alilen = len(row["qseq"])
        qfrag, sfrag = [], []
        for i in range(alilen):
            current_region_type = region_type(q_cur, start_codon,stop_codon)
            if current_region_type!= last_region_type:
                partitioned_alis.append([row["qseqid"], row["sseqid"], q_fr_start, q_cur-1, s_fr_start, s_cur-1, qframe,row["qlen"], "".join(qfrag), "".join(sfrag),last_region_type])
                qfrag, sfrag = [], []
                q_fr_start = q_cur
                s_fr_start = s_cur
            qfrag.append(row["qseq"][i])
            sfrag.append(row["sseq"][i])
            if row["qseq"][i] != "-":
                q_cur += 3
            if row["sseq"][i] != "-":
                s_cur += 1
            last_region_type = current_region_type
        partitioned_alis.append([row["qseqid"], row["sseqid"], q_fr_start, q_cur-1, s_fr_start, s_cur-1, qframe,row["qlen"], "".join(qfrag), "".join(sfrag),last_region_type])
        type_broken.extend(partitioned_alis)
    type_broken_df = pd.DataFrame(type_broken, columns =["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "qlen", "qseq", "sseq", "region_type"])
    return type_broken_df #, columns = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "qseq", "sseq"]