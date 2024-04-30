def total_area_covered(df, start_col= "bounded_dom_start", end_col= "bounded_dom_end"):
    cov = set()
    for index,row in df.iterrows():
        cov = cov.union(set(range(row[start_col],row[end_col]+1)))
    return len(cov)

def calc_bounded_dom_cov(alidf, domdata, col4domcov, colstart, colend):
    """This function would be used for finding the domain coverage in a specific region
    that is specified by colstart and colend"""
    ali_dom = alidf.merge(domdata, left_on = col4domcov, right_on = "protein_id", how = "left").drop(columns= ["protein_id"])
    #The next two lines are for handling the sequences without a domain 
    ali_dom["domstart"] = ali_dom["domstart"].fillna(10).astype(int)
    ali_dom["domend"] = ali_dom["domend"].fillna(0).astype(int)
    ali_dom["bounded_dom_start"] = ali_dom[["domstart", colstart]].max(axis=1)
    ali_dom["bounded_dom_end"] = ali_dom[["domend", colend]].min(axis=1)
    ali_with_cov = ali_dom.groupby(list(alidf.columns)).apply(total_area_covered).reset_index()
    return ali_with_cov