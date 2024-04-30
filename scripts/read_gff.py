from gff_cols import gff_cols
import pandas as pd

def read_gff(path):
    """This function gets the path to the gff file as input and returns a dataframe 
    by reading the gff file. It also names the columns of the gff file. """
    gff = pd.read_csv(path, sep="\t", header=None, comment="#", names = gff_cols )
    return gff

