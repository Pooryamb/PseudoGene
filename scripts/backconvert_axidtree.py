#!/usr/bin/python3
import os
import glob
import json

def map_axidtree2geneid(axid_tree, id_mapping):
    non_alnums = ['(', ')', ',', ':']
    for char in non_alnums:
        axid_tree = axid_tree.replace(char, f"\n{char}\n")
    parts = axid_tree.split("\n")
    geneid_parts = []
    for part in parts:
        if part.startswith("A"):
            geneid_parts.append(id_mapping[part.strip()])
            if part.strip()=="A0":
                geneid_parts.append("{test}")
        else:
            geneid_parts.append(part)
    return "".join(geneid_parts) + "\n"

phyml_trees = glob.glob("../data/manual_checking/*/*/seq4dnds/orths.PHYLIP_phyml_tree.txt")
for tree in phyml_trees:
    json_file_path = tree.replace(".PHYLIP_phyml_tree.txt", "_idmap.json")

    nwk_tree = open(tree).read()
    id_mapping = json.loads(open(json_file_path).read())

    new_geneid_path = tree.replace(".PHYLIP_phyml_tree.txt", "_geneid.nwk")

    with open(new_geneid_path, 'w') as outfile:
        outfile.write(map_axidtree2geneid(nwk_tree, id_mapping))
