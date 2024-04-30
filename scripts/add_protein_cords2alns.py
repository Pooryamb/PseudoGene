#!/usr/bin/python3
import pandas as pd
import numpy as np

def add_protein_cords2alns(blastx_aln):
    df_aln = blastx_aln.copy()
    pstart, pend, protlen = "qstart_p", "qend_p", "protlen_q"
    """The input must have the start_codon and stop_codon"""
    #df_aln = blastx_aln.merge(start_stop_pos, left_on = "qseqid", right_on = "gene_id").drop(columns= ["gene_id"])
    df_aln[pstart] = np.ceil((df_aln["qstart"] - df_aln["start_codon"] + 1 )/3).astype(int)
    df_aln[pend  ] = np.ceil((df_aln["qend"] - df_aln["start_codon"] + 1 )/3).astype(int)
    df_aln[protlen] = (df_aln["stop_codon"] - df_aln["start_codon"] + 1)//3
    return df_aln
