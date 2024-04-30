#!/usr/bin/python3
import glob
import os
import json
from fasta2dict import fasta2dict


nal_files_path = glob.glob("../data/manual_checking/*/*/seq4dnds/*.nal")

for file_path in nal_files_path:
    seq_dict = fasta2dict(file_path)
    new2tritryp_mapping = {f"A{i}":list(seq_dict.keys())[i] for i in range(len(seq_dict))}
    tritryp2new_mapping = {y:x for (x,y) in new2tritryp_mapping.items()}
    
    num_of_seqs = len(seq_dict)
    seq_len = len(list(seq_dict.values())[0])
    
    with open(file_path.replace(".nal", ".PHYLIP"), 'w') as phylip_file, \
         open(file_path.replace(".nal", "_idmap.json"), 'w') as phylip_id_mapping:
        json.dump(new2tritryp_mapping, phylip_id_mapping)
        
        phylip_file.write(f" {num_of_seqs} {seq_len}\n")
        for i, seq in enumerate(seq_dict.values()):
            new_seq_id = f"A{i}"
            phylip_file.write(f"{new_seq_id.ljust(10)}{seq}\n")
