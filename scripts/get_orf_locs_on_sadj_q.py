import pandas as pd

def get_orf_locs_on_sadj_q(row):
    qseq, qstart, qend, start_codon, stop_codon = row["qseq"], row["qstart"], row["qend"], row["start_codon"], row["stop_codon"]
    coding_start = 0
    coding_end = 0
    includes_start_codon, includes_stop_codon, found_start_codon, found_stop_codon= True, True, False, False
    if qstart > start_codon:
        includes_start_codon = False
        found_start_codon = False
    if qend < stop_codon:
        includes_stop_codon = False
        found_stop_codon = False

    q_loc_nuc = qstart - 1
    q_loc_pr = 0

    for letter in qseq:
        if letter == "-":
            pass
        elif letter == "/":
            q_loc_nuc -= 1
        elif letter == "\\":
            q_loc_nuc += 1
        else:
            q_loc_pr +=1
            q_loc_nuc += 3
        if (includes_start_codon and not(found_start_codon)) and (q_loc_nuc >= start_codon):
            coding_start = q_loc_pr
            found_start_codon = True
        if (includes_stop_codon and not(found_stop_codon)) and (q_loc_nuc >= stop_codon):
            coding_end = q_loc_pr
            found_stop_codon = True
            break
    return pd.Series([coding_start, coding_end])