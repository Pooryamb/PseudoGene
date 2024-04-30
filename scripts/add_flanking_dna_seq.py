#!/usr/bin/python3
from fasta2dict import fasta2dict
from select_first_name_in_header import select_first_name_in_header

selected_frag_dict = fasta2dict("../data/cds_and_flanking.fasta", select_first_name_in_header)


def add_5prime_end(row, frame_type):
    if row["start_codon_on_sadj_pr"]!=0:
        return ""
    diff = row["qstart"] - row["start_codon"]
    diff_div_3 = diff - diff%3
    wholeseq = selected_frag_dict[row["qseqid"]]
    if frame_type=="ori":
        start_seq = row["start_codon"]-1
        if diff%3==1:
            end_seq = row["start_codon"]+ diff_div_3 -1
        elif diff%3==2:
            end_seq = row["start_codon"]+ diff_div_3 -1 +3
        else:
            end_seq = row["qstart"] - 1
            #print("difference is not valid")
    elif frame_type=="new":
        start_seq =row["qstart"]-diff_div_3-1
        end_seq = row["qstart"]-1
    return wholeseq[start_seq: end_seq]

def add_3prime_end(row, frame_type):
    if row["stop_codon_on_sadj_pr"]!=0:
        return ""
    diff = row["stop_codon"] - row["qend"]
    diff_div_3 = diff - diff%3
    wholeseq = selected_frag_dict[row["qseqid"]]
    if frame_type=="ori":
        end_seq = row["stop_codon"]
        if diff%3==1:
            start_seq = end_seq - diff_div_3 
        elif diff%3==2:
            start_seq = end_seq - diff_div_3 -3
        else:
            start_seq =row["qend"]
    elif frame_type=="new":
        start_seq =row["qend"]
        end_seq = start_seq + diff_div_3
    return wholeseq[start_seq: end_seq]
