#!/usr/bin/python3

import json
import glob
import os
import sys
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from read_ipr import read_ipr
from calc_dom_cov import total_area_covered

dom_trunc_per_inc_thresh = 0.2
dom_trunc_len_inc_thresh = 20
coding_dom_per_inc_thresh = 0.2
coding_dom_len_inc_thresh = 20

relax_json_paths = glob.glob("../data/manual_checking/*/*/orths.nal.RELAX.json")
relax_analysis = []
for json_path in relax_json_paths:
    with open(json_path) as json_file:
        gene_id = os.path.basename(os.path.dirname(json_path))
        relax_dict = json.loads(json_file.read())
        p_value, k = relax_dict["test results"]["p-value"], relax_dict["test results"]["relaxation or intensification parameter"]
        relax_analysis.append([gene_id, k, p_value])
relax_analysis_df = pd.DataFrame(relax_analysis, columns = ["pseudogene", "k(RELAX)", "p_value(RELAX)"])

all_pgs_paths = glob.glob("../data/manual_checking/*/*/")
all_pf_pgs_paths = glob.glob("../data/manual_checking/pf/*/")

# This part is for finding the difference between the number of rbhs of pseudogene and their templates
# orthos_tbl_path = sys.argv[1]
orthos_tbl_path = "../data/rbh_tb927_trypanosomes.tsv"
rbh_df = pd.read_csv(orthos_tbl_path, sep="\t")
source_org, closest_org = rbh_df.columns[:2]
pg_subj_orths_count = []
for pg_path in all_pgs_paths:
    pg, subj, _, pident = open(f"{pg_path}/nonbroken_ali.tsv").readlines()[1].split()[:4]
    if pg in list(rbh_df[source_org]):
        pg_orth_count = (~rbh_df[rbh_df[source_org] == pg].isna()).iloc[0].sum() - 1
        pg_in_closeorg = not(rbh_df[rbh_df[source_org] == pg][closest_org].isna().iloc[0])
    else:
        pg_orth_count = 0
        pg_in_closeorg = False
    
    if subj in list(rbh_df[source_org]):
        subj_orth_count = (~rbh_df[rbh_df[source_org] == subj].isna()).iloc[0].sum() - 1
    else:
        subj_orth_count = 0
    pg_subj_orths_count.append([pg, subj,pident,  pg_orth_count, subj_orth_count, pg_in_closeorg])
pg_subj_orths_count_df = pd.DataFrame(pg_subj_orths_count, columns = ["pseudogene", "subject", "pident", "pg_orth_count", "subj_orth_count", "pg_isin_close_organism"])

# The following part is for evaluating changes in domain coverage that corresponds to coding regions or truncates a domain

original_domains = read_ipr("../data/ipr_proteome.tsv")
original_domcov = original_domains.groupby("protein_id").apply(lambda grp: total_area_covered(grp, "domstart", "domend")).reset_index().rename(columns = {0: "dom_cov"})


truncated_domains_all_pgs = []
contains_inc_in_coding_cov = []
for pg_path in all_pf_pgs_paths:
    corrected_doms = pd.read_csv(f"{pg_path}/correctedseq_doms.tsv", sep="\t")
    pg_id = os.path.basename(os.path.normpath(pg_path))
    corrected_doms["original_dom_length"] = np.maximum(corrected_doms["domend"] - corrected_doms["domstart"] + 1, 0)
    corrected_doms["start_codon_trun_dom_start"] = corrected_doms[["start_codon_pos", "domstart"]].max(axis=1)
    corrected_doms["stop_codon_trun_dom_end"] = corrected_doms[["stop_codon_pos", "domend"]].min(axis=1)
    start_codon_pos, stop_codon_pos = corrected_doms.iloc[0]["start_codon_pos"], corrected_doms.iloc[0]["stop_codon_pos"]
    corrected_doms["truncated_dom_length"] = np.maximum(corrected_doms["stop_codon_trun_dom_end"] - corrected_doms["start_codon_trun_dom_start"] + 1, 0)
    truncated_doms = corrected_doms[((corrected_doms["truncated_dom_length"]/corrected_doms["original_dom_length"] <= (1 - dom_trunc_per_inc_thresh)) &
                   (corrected_doms["truncated_dom_length"]/corrected_doms["original_dom_length"] >= dom_trunc_per_inc_thresh)) &
                   (corrected_doms["original_dom_length"] - corrected_doms["truncated_dom_length"] >= dom_trunc_len_inc_thresh )
                  ]

    if truncated_doms.shape[0] != 0 :
        truncated_doms["dom_ident"] = truncated_doms["db_acc"] + "_" + corrected_doms["domstart"].astype(str)  + "_" + corrected_doms["domend"].astype(str)
        truncated_doms_summary = ",".join(truncated_doms["dom_ident"])
        truncated_domains_all_pgs.append([pg_id, True, truncated_doms_summary])
    else:
        truncated_domains_all_pgs.append([pg_id, False, ""])
    
    ori_dom_cov_on_pg_df = original_domcov[original_domcov["protein_id"] == pg_id]
    if ori_dom_cov_on_pg_df.empty:
        ori_dom_cov_on_pg = 0
    else:
        ori_dom_cov_on_pg = ori_dom_cov_on_pg_df.iloc[0]["dom_cov"]
                                        
    coding_dom_cov = total_area_covered(corrected_doms, "start_codon_trun_dom_start", "stop_codon_trun_dom_end")
    
    if (coding_dom_cov - ori_dom_cov_on_pg >= coding_dom_len_inc_thresh) and (ori_dom_cov_on_pg == 0 or (coding_dom_cov/ori_dom_cov_on_pg >= (1 + coding_dom_per_inc_thresh))):
        contains_inc_in_coding_cov.append([pg_id, True, coding_dom_cov, ori_dom_cov_on_pg,start_codon_pos, stop_codon_pos ])
    else:
        contains_inc_in_coding_cov.append([pg_id, False, coding_dom_cov, ori_dom_cov_on_pg,start_codon_pos, stop_codon_pos])

truncated_domains_all_pgs = pd.DataFrame(truncated_domains_all_pgs, columns = ["pseudogene", "contains_truncated_domains (only pf)", "truncated_domains (only pf)"])
contains_inc_in_coding_cov = pd.DataFrame(contains_inc_in_coding_cov, columns = ["pseudogene", "sign_inc_in_coding_domcov", "domcov_in_coding (only pf)", "domcov_in_original (only pf)", "start_codon_pos", "stop_codon_pos"])


all_pgs_df = pd.DataFrame({"pseudogene": [os.path.basename(os.path.normpath(x)) for x in all_pgs_paths]})
all_pgs_summary = all_pgs_df.merge(pg_subj_orths_count_df, how = "left", on = "pseudogene").merge(relax_analysis_df, how = "left", on = "pseudogene").merge(truncated_domains_all_pgs, how = "left", on = "pseudogene").merge(contains_inc_in_coding_cov, how = "left", on = "pseudogene") 
nf_pgs = [os.path.basename(os.path.normpath(x)) for x in glob.glob("../data/manual_checking/nf/*/")]
all_pgs_summary["positive frame"] = all_pgs_summary.apply(lambda row: row["pseudogene"] not in nf_pgs, axis= 1)
all_pgs_summary.to_csv("../data/summary_of_pgs.tsv", sep = "\t", index = None)
