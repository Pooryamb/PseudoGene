#!/usr/bin/python3
def break_sadj_qname_pf(df):
    dfc = df.copy()
    dfc[["___", "subject", "start_info", "end_info"]] = dfc["protein_id"].str.split("__", expand = True)
    dfc[["start_frame", "start_codon_pos"]] = dfc["start_info"].str.split("_", expand= True).astype(int)
    dfc[["end_frame", "stop_codon_pos"]] = dfc["end_info"].str.split("_", expand= True).astype(int)
    all_input_cols = list(df.columns)
    all_input_cols.remove("protein_id")
    all_input_cols.remove("query")
    selected_cols = ["query", "subject", "start_frame", "start_codon_pos", "end_frame", "stop_codon_pos"] + all_input_cols
    return dfc[selected_cols]