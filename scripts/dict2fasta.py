#!/usr/bin/python3
def dict2fasta(seqdict, outputpath):
    fasta_file = open(outputpath, 'w')
    for seqid, seq in seqdict.items():
        fasta_file.write(f">{seqid}\n{seq}\n")
    fasta_file.close()
    