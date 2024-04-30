import pandas as pd

def has_fs(aln):
    """This function returns a boolean to show if the lines have frameshifts"""
    has_fs = (aln['qseq'].str.contains(r'/|\\', regex=True)) | (aln["qframe"]!=1)
    return has_fs
    
def is_trun(aln):
    """This function returns a boolean to show if the lines are truncated"""
    qstart = aln[['qstart', 'qend']].min(axis=1)
    qend = aln[['qstart', 'qend']].max(axis=1)
    start_codon = aln[['start_codon', 'stop_codon']].min(axis=1)
    stop_codon = aln[['start_codon', 'stop_codon']].max(axis=1)
    is_trun = (qstart<start_codon) | (qend>stop_codon)
    return is_trun 

def select_fs_or_trun_alis(aln):
    """This function assumes all the rows cover at least a part of the coding part of the query."""
    trun_lines = is_trun(aln)
    fs_lines = has_fs(aln)
    
    selected = aln[trun_lines|fs_lines]
    return selected.reset_index(drop=True)
    