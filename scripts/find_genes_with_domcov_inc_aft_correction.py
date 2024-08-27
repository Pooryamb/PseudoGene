#!/usr/bin/python3

from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from read_ipr import read_ipr
from calc_dom_cov import total_area_covered, calc_bounded_dom_cov
from has_sig_domcov_inc import has_sig_domcov_inc
from break_sadj_qname_pf import break_sadj_qname_pf
from break_dmndali2frames_and_types import break_dmndali2frames_and_types
from break_dmndali2frames import break_dmndali2frames
from reverse_frame import reverse_frame



ipr_path = "../data/ipr_proteome.tsv"
ipr_ori_prot = read_ipr(ipr_path)

dmnd_alis = pd.read_csv("../data/genome_vs_proteome_fmt6_alis.txt", sep="\t")

ori_proteome_dom_cov = ipr_ori_prot.groupby("protein_id").apply(lambda x: total_area_covered(x, "domstart", "domend")).reset_index().rename(columns={0:"dom_cov_on_original_protein"})

ipr_sadj_prot = read_ipr("../data/ipr_sadj_seq.tsv")

# this cell is for reading the location of start/stop codon on each query fragment
start_stop_pos_ori = pd.read_csv("../data/start_stop_pos_ori.tsv", sep="\t")

#First, we sort the subject adjusted sequences based on their increase in the domain coverage and consider the gene with the highest domain coverage as the main corrected version.
ipr_sadj_prot["query"] = ipr_sadj_prot["protein_id"].str.split("__", expand =True)[0]
corrected_proteome_dom_cov = ipr_sadj_prot.groupby(["protein_id", "query"]).apply(lambda x: total_area_covered(x, "domstart", "domend")).reset_index().rename(columns={0:"dom_cov_on_corrected_protein"})
corrections_with_max_domcov = corrected_proteome_dom_cov.sort_values(by=["protein_id", "dom_cov_on_corrected_protein"], ascending = [True, False]).drop_duplicates("query").reset_index(drop=True)

#Here, we want to select the proteins whose domain coverage increases significantly after correction
original_vs_corrected_domcov = corrections_with_max_domcov.merge(ori_proteome_dom_cov, left_on = "query", right_on = "protein_id", how = "left").drop(columns = ["protein_id_y"]).rename(columns={"protein_id_x": "protein_id"})
original_vs_corrected_domcov["dom_cov_on_original_protein"] = original_vs_corrected_domcov["dom_cov_on_original_protein"].fillna(0).astype(int)
dom_cov_len_inc_thresh = 20
dom_cov_per_inc_thresh = 0.2
sig_inc_dom_cov_after_correction = has_sig_domcov_inc(original_vs_corrected_domcov, dom_cov_len_inc_thresh, dom_cov_per_inc_thresh,  "dom_cov_on_original_protein", "dom_cov_on_corrected_protein")

ipr_sadj_prot = sig_inc_dom_cov_after_correction[["protein_id"]].merge(ipr_sadj_prot)

#First, we take care of query-subject alignments that query aligns with subject with a positive frame
ipr_sadj_prot_p = ipr_sadj_prot[~(ipr_sadj_prot["protein_id"].str.contains("__negf"))]
if not(ipr_sadj_prot_p.empty):
    ipr_sadj_prot_p = break_sadj_qname_pf(ipr_sadj_prot_p)

    dmnd_alis_pg_p = dmnd_alis.merge(ipr_sadj_prot_p[["query", "subject"]].drop_duplicates(), left_on =["qseqid", "sseqid"] , right_on =["query", "subject"])

    broken_alis_pg_p = break_dmndali2frames_and_types(dmnd_alis_pg_p, start_stop_pos_ori)

    for gene_id in set(ipr_sadj_prot_p["query"]):
        dir_path = f"../data/manual_checking/pf/{gene_id}"
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        gene_dmnd_ali = dmnd_alis_pg_p[dmnd_alis_pg_p["qseqid"] == gene_id]
        gene_dmnd_ali.to_csv(f"{dir_path}/nonbroken_ali.tsv", sep="\t", index=None)
        subj_id = gene_dmnd_ali.iloc[0]["sseqid"]
        gene_broken_dmnd_ali = broken_alis_pg_p[broken_alis_pg_p["qseqid"] == gene_id]
        gene_broken_dmnd_ali.to_csv(f"{dir_path}/broken_ali.tsv", sep="\t", index=None)
        gene_ipr_doms = ipr_ori_prot[ipr_ori_prot["protein_id"] == gene_id]
        gene_ipr_doms.to_csv(f"{dir_path}/ori_f_doms.tsv", sep="\t", index=None)
        sub_gene_ipr_doms = ipr_ori_prot[ipr_ori_prot["protein_id"] == subj_id]
        sub_gene_ipr_doms.to_csv(f"{dir_path}/subject_doms.tsv", sep="\t", index=None)
        corrected_ipr_doms = ipr_sadj_prot_p[ipr_sadj_prot_p["query"] == gene_id]
        corrected_ipr_doms.to_csv(f"{dir_path}/correctedseq_doms.tsv", sep="\t", index=None)
    
#Now, we take care of query-subject alignments that query aligns with subject with a negative frame
#For those aligning with a negative frame, having the broken alignments is not helpful
ipr_sadj_prot_n = ipr_sadj_prot[ipr_sadj_prot["protein_id"].str.contains("__negf")]
if not(ipr_sadj_prot_n.empty):
    ipr_sadj_prot_n["subject"] = ipr_sadj_prot["protein_id"].str.split("__", expand=True)[1]
    ipr_sadj_prot_n = ipr_sadj_prot_n[['query', 'subject', 'db', 'db_acc', 'db_desc', 'ipr_desc', 'domstart', 'domend']]

    dmnd_alis_pg_n = dmnd_alis.merge(ipr_sadj_prot_n[["query", "subject"]].drop_duplicates(), left_on =["qseqid", "sseqid"] , right_on =["query", "subject"]).drop(columns = ["query", "subject"])
    reversed_aln = reverse_frame(dmnd_alis_pg_n, ["qlen", ["qstart", "qend"]],{"qframe":{x:-x for x in [1,2,3,-1,-2,-3]}} )
    broken_rev_alis_pg_n = break_dmndali2frames(reversed_aln)

    for gene_id in set(ipr_sadj_prot_n["query"]):
        dir_path = f"../data/manual_checking/nf/{gene_id}"
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        gene_dmnd_ali = dmnd_alis_pg_n[dmnd_alis_pg_n["qseqid"] == gene_id]
        gene_dmnd_ali.to_csv(f"{dir_path}/nonbroken_ali.tsv", sep="\t", index=None)
        gene_broken_dmnd_ali = broken_rev_alis_pg_n[broken_rev_alis_pg_n["qseqid"] == gene_id]
        gene_broken_dmnd_ali.to_csv(f"{dir_path}/broken_rev_ali.tsv", sep="\t", index=None)
        subj_id = gene_dmnd_ali.iloc[0]["sseqid"]
        
        gene_ipr_doms = ipr_ori_prot[ipr_ori_prot["protein_id"] == gene_id]
        gene_ipr_doms.to_csv(f"{dir_path}/ori_f_doms.tsv", sep="\t", index=None)
        sub_gene_ipr_doms = ipr_ori_prot[ipr_ori_prot["protein_id"] == subj_id]
        sub_gene_ipr_doms.to_csv(f"{dir_path}/subject_doms.tsv", sep="\t", index=None)
        corrected_ipr_doms = ipr_sadj_prot_n[ipr_sadj_prot_n["query"] == gene_id]
        corrected_ipr_doms.to_csv(f"{dir_path}/correctedseq_doms.tsv", sep="\t", index=None)
