def has_sig_domcov_inc(ali, len_thresh, perc_thresh, col1domcov = "q_dom_cov", col2domcov = "s_dom_cov"):
    dom_cov_diff = ali[col2domcov] - ali[col1domcov]
    cond1 = (dom_cov_diff >= len_thresh)
    cond2 = (dom_cov_diff/ali[col1domcov] >= perc_thresh )
    return ali[cond1 & cond2]