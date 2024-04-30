def find_domain_contatining_genes_trun_by_start_stop_frame(df, cond1cols, cond2cols, lenthresh, percthresh):
    domlen_fl = df[cond1cols[1]] -df[cond1cols[0]] +1
    domlen_trun = df[cond2cols[1]] -df[cond2cols[0]] +1
    df_trun = df[((domlen_fl - domlen_trun) >= lenthresh) &(domlen_fl/domlen_trun >= 1+ percthresh) ]
    return df_trun