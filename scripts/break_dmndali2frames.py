#!/usr/bin/python3
import pandas as pd
def break_dmndali2frames(alidf):
    """Note: This function works only on alignments with a positive frame
    it takes a row as input and returns a list of rows as output"""
    
    fs_dict = {1:3, 2:1, 3:2} #This is frame change after seeing a  forward slash
    bs_dict = {1:2, 2:3, 3:1} #This is frame change after seeing a backward slash
    frame_broken = []
   
    for index,row in alidf.iterrows():
        partitioned_alis = []   #Each section will be represented by "qstart", "qend", "sstart", "send", "qframe", "qseq", "sseq"
        q_newfr_start,s_newfr_start,qframe = row["qstart"], row["sstart"], row["qframe"]
        q_cur, s_cur = q_newfr_start, s_newfr_start
        #q_cur and s_cur is for keeping the record of the position of the cursor on query and subject
        alilen = len(row["qseq"])
        qfrag, sfrag = [], []
        for i in range(alilen):
            if (row["qseq"][i] == "/") or (row["qseq"][i] == "\\"):
                partitioned_alis.append([row["qseqid"], row["sseqid"], q_newfr_start, q_cur-1, s_newfr_start, s_cur-1, qframe, row["qlen"], "".join(qfrag), "".join(sfrag)])
                qfrag, sfrag = [], []
                if (row["qseq"][i] == "/"):
                    qframe = fs_dict[qframe]
                    q_cur -= 1
                else:
                    qframe = bs_dict[qframe]
                    q_cur += 1
                q_newfr_start = q_cur
                s_newfr_start = s_cur
            else:
                qfrag.append(row["qseq"][i])
                sfrag.append(row["sseq"][i])
                if row["qseq"][i] != "-":
                    q_cur += 3
                if row["sseq"][i] != "-":
                    s_cur += 1
        partitioned_alis.append([row["qseqid"], row["sseqid"], q_newfr_start, q_cur-1, s_newfr_start, s_cur-1, qframe,row["qlen"], "".join(qfrag), "".join(sfrag)])
        frame_broken.extend(partitioned_alis)
    frame_broken_df = pd.DataFrame(frame_broken, columns =["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "qlen", "qseq", "sseq"] )
    return frame_broken_df #, columns = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe","qlen", "qseq", "sseq"]
    
    