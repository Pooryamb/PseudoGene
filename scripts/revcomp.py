def revcomp(seq):
    seqcomp = seq.replace("A","t").replace("T","a").replace("C","g").replace("G","c").upper()
    return seqcomp[::-1]