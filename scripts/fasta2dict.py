def fasta2dict(fastapath, header_processor = lambda x: x):
    seq_dict = {}
    parts = open(fastapath).read().strip().lstrip(">").split("\n>")
    for part in parts:
        lines = part.split("\n")
        seq_id = header_processor(lines[0])
        seq = "".join(lines[1:])
        seq_dict[seq_id] = seq
    return seq_dict