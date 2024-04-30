def get_poses_reversed_frame(locationdf):
    locdf = locationdf.copy()
    locdf.loc[locdf["strand"]=="+", "start_codon"] = locdf["frag_end"] - locdf["cds_end"] + 1
    locdf.loc[locdf["strand"]=="+", "stop_codon"] = locdf["frag_end"] - locdf["cds_start"] + 1

    locdf.loc[locdf["strand"]=="-", "start_codon"] = locdf["cds_start"] - locdf["frag_start"] + 1
    locdf.loc[locdf["strand"]=="-", "stop_codon"] = locdf["cds_end"] - locdf["frag_start"] + 1

    locdf["start_codon"] = locdf["start_codon"].astype(int)
    locdf["stop_codon"] = locdf["stop_codon"].astype(int)

    start_stop_codon = locdf[["gene_id", "start_codon", "stop_codon"]]
    return start_stop_codon


def get_poses_original_frame(locationdf):
    locdf = locationdf.copy()
    locdf.loc[locdf["strand"]=="+", "start_codon"] = locdf["cds_start"] - locdf["frag_start"] + 1
    locdf.loc[locdf["strand"]=="+", "stop_codon"] = locdf["cds_end"] - locdf["frag_start"] + 1

    locdf.loc[locdf["strand"]=="-", "start_codon"] = locdf["frag_end"] - locdf["cds_end"] + 1
    locdf.loc[locdf["strand"]=="-", "stop_codon"] = locdf["frag_end"] - locdf["cds_start"] + 1

    locdf["start_codon"] = locdf["start_codon"].astype(int)
    locdf["stop_codon"] = locdf["stop_codon"].astype(int)

    start_stop_codon = locdf[["gene_id", "start_codon", "stop_codon"]]
    return start_stop_codon
