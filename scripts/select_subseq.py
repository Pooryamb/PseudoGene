from revcomp import revcomp

def select_subseq(seqdict, locdf):
    """This function takes a seqdict as input, and a dataframe containing 
    the gene name, contig name, start, and end on fragment and it returns 
    a seqdict of the selected regions"""
    outputseqdict = {}
    for index,row in locdf.iterrows():
        seq = seqdict[row["contig"]][row["frag_start"]-1:row["frag_end"]]
        if row["strand"]=="-":
            seq = revcomp(seq)
        outputseqdict[row["gene_id"]] = seq
    return outputseqdict