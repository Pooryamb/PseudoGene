#NOTE: This script works only on positive frames.

import pandas as pd
from break_dmndali2frames import break_dmndali2frames
from break_ali2genomictypes import break_ali2genomictypes

def break_dmndali2frames_and_types(alidf, start_stop_locations):
    
    frame_broken_df = break_dmndali2frames(alidf)
    frame_broken_df = frame_broken_df.merge(start_stop_locations, left_on = "qseqid", right_on = "gene_id").drop(columns= ["gene_id"])
    
    type_broken_df = break_ali2genomictypes(frame_broken_df)
    return type_broken_df