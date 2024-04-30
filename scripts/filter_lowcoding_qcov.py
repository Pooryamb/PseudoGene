import pandas as pd
from overlap import overlap

def filter_lowcoding_qcov(ali_df, qcoding_cov = 0.5):
    """Each line of the input table should also contain the 
    start_codon and stop_codon locations"""
    qstart = ali_df[["qstart", "qend"]].min(axis=1)
    qend = ali_df[["qstart", "qend"]].max(axis=1)
    df_temp = pd.DataFrame({"qstart": qstart, "qend": qend, 
                       "start_codon":ali_df["start_codon"], 
                       "stop_codon":ali_df["stop_codon"]})
    coding_overlap = overlap(df_temp, ["qstart", "qend"],["start_codon", "stop_codon"])
    coding_overlap = coding_overlap/(df_temp["stop_codon"] - df_temp["start_codon"]+ 1)
    ali_df["qcoding_cov"] = coding_overlap
    filtered = ali_df[ali_df["qcoding_cov"] >= qcoding_cov]
    return filtered