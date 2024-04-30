#!/usr/bin/python3

import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from filter_lowcoding_qcov import filter_lowcoding_qcov
from select_fs_or_trun_alis import select_fs_or_trun_alis
from read_ipr import read_ipr
from add_protein_cords2alns import add_protein_cords2alns
from calc_dom_cov import calc_bounded_dom_cov, total_area_covered
from has_sig_domcov_inc import has_sig_domcov_inc
from find_end_frame import find_end_frame
from get_orf_locs_on_sadj_q import get_orf_locs_on_sadj_q
from format_ali_seq import format_ali_seq
from get_flanking_prseq_and_info import get_5_prseq_and_info, get_3_prseq_and_info
from dict2fasta import dict2fasta

def sort_two_cols(df, col1, col2):
    df_copy = df.copy()
    col1_sorted = df[[col1,col2]].min(axis=1)
    col2_sorted = df[[col1,col2]].max(axis=1)
    df_copy[col1], df_copy[col2] = col1_sorted, col2_sorted
    return df_copy

start_stop_loc = pd.read_csv("../data/start_stop_pos_ori.tsv", sep="\t")

ipr_dom = read_ipr("../data/ipr_proteome.tsv")

alignments = pd.read_csv(f"../data/genome_vs_proteome_fmt6_alis.txt", sep = "\t")
alignments = alignments.merge(start_stop_loc, left_on = "qseqid", right_on = "gene_id").drop(columns=["gene_id"])
alignments = filter_lowcoding_qcov(alignments)
alignments = alignments[alignments["qseqid"] != alignments["sseqid"]]
alignments = select_fs_or_trun_alis(alignments)
alignments = add_protein_cords2alns(alignments)
alignments = sort_two_cols(alignments,"qstart_p", "qend_p" )
alignments = calc_bounded_dom_cov(alignments, ipr_dom, "sseqid", "sstart","send").rename(columns={0:"s_dom_cov"})
alignments = alignments.drop_duplicates(["qseqid", "qseq"])

proteome_dom_cov = ipr_dom.groupby("protein_id").apply(lambda x: total_area_covered(x, "domstart", "domend")).reset_index().rename(columns={0:"dom_cov_on_protein"})

#_p and _n denote positive and negative frames, respectively.
alignments_p = alignments[alignments["qframe"]>0]
alignments_p = calc_bounded_dom_cov(alignments_p, ipr_dom, "qseqid", "qstart_p","qend_p").rename(columns={0:"q_dom_cov"})

alignments_n = alignments[alignments["qframe"]<0]
alignments_n = alignments_n.merge(proteome_dom_cov, left_on = "qseqid", right_on = "protein_id", how= "left").drop(columns=["protein_id"])
alignments_n = alignments_n.rename(columns = {"dom_cov_on_protein":"q_dom_cov"})
alignments_n["q_dom_cov"] = alignments_n["q_dom_cov"].fillna(0).astype(int)

#These numbers specify the threshold for domain coverage increase to investigate if it is a pseudogene.

dom_cov_len_inc_thresh = 20
dom_cov_per_inc_thresh = 0.2

sig_alignments_n = has_sig_domcov_inc(alignments_n, dom_cov_len_inc_thresh, dom_cov_per_inc_thresh)

sig_alignments_p = has_sig_domcov_inc(alignments_p, dom_cov_len_inc_thresh, dom_cov_per_inc_thresh)
# For psuedogenes predicted on the same direction as the current gene (with a positive frame), we need to know the frame
# at the start and the end of the alignment. But for the negative frames, we don't need such a thing because it is not possible
# to have a protein sequence from negative and positive frames.

sig_alignments_p["qframe_end"] = sig_alignments_p.apply(lambda x: find_end_frame(x["qseq"], x["qframe"]), axis=1)
sig_alignments_p[["start_codon_on_sadj_pr", "stop_codon_on_sadj_pr"]] = sig_alignments_p.apply(get_orf_locs_on_sadj_q, axis=1)

adjusted_seq_dict = {}

for index,row in sig_alignments_p.iterrows():
    infos_3, seqs_3 = get_3_prseq_and_info(row)
    infos_5, seqs_5 = get_5_prseq_and_info(row)
    for i in range(len(infos_3)):
        for j in range(len(infos_5)):
            fr_3, end_3  = infos_3[i]
            fr_5, start_5 = infos_5[j]
            seq_3 = seqs_3[i]
            seq_5 = seqs_5[j]
            adjusted_seq = seq_5 + row["qseq"] + seq_3
            start_loc = start_5
            if seq_3 != "":
                end_loc = len(seq_5) + end_3
            else:
                end_loc = len(seq_5) + len(format_ali_seq(row["qseq"]))+ end_3
            # The header of each sequence contains:
            # qseqid, sseqid, 5' frame, original start codon loc, 3' frame, original stop codon loc
            seq_name = row["qseqid"] + "__" +  row["sseqid"] + f"__{fr_5}_{start_5}__{fr_3}_{end_loc}"
            seq = format_ali_seq(seq_5 + row["qseq"] + seq_3)
            adjusted_seq_dict[seq_name] = seq

for index,row in sig_alignments_n.iterrows():
    qseqid = row["qseqid"]
    sseqid = row["sseqid"]
    seq_name = row["qseqid"] + "__" +  row["sseqid"] + "__negf"
    seq = format_ali_seq(row["qseq"])
    adjusted_seq_dict[seq_name] = seq

dict2fasta(adjusted_seq_dict, "../data/subject_adjusted_proteinseq.fasta")