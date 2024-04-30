import pandas as pd

ipr_header = ['protein_id', 'seq_md5', 'protlen', 'db', 'db_acc', 'db_desc', 
    'domstart', 'domend', 'score', 'ismatch', 'data', 'ipr_acc', 'ipr_desc']

def read_ipr(path):
    df = pd.read_csv(path, sep="\t", header=None, names= ipr_header)
    selected_dbs = ["Gene3D", "Pfam", "SUPERFAMILY","SMART","PRINTS", "ProSiteProfiles", "FunFam", "CDD", "ProSitePatterns", "PIRSF", "Hamap", "SFLD"]
    df = df[df["db"].isin(selected_dbs)]
    selected_cols= ['protein_id', 'db', 'db_acc', 'db_desc', 'ipr_desc', 'domstart', 'domend']
    return df[selected_cols].reset_index(drop=True)
    