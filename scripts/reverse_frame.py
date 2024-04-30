import pandas as pd

def reverse_frame(alndf, coords2change, column_change_dict):
    """The second argument should be in the format of a list where its first element shows the column containing the 
    length of the query and its second element should be a list of start and ends """
    len_col, coords_cols = coords2change
    rev_df = alndf.copy()
    for column, mapping_dict in column_change_dict.items():
        rev_df[column] = alndf[column].map(mapping_dict)
    for i in range(0,len(coords_cols),2):
        start_col, end_col = coords_cols[i], coords_cols[i+1]
        rev_df[start_col] = alndf[len_col] - alndf[start_col] + 1
        rev_df[end_col]   = alndf[len_col] - alndf[end_col]   + 1

    return rev_df