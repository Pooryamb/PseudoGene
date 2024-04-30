from add_flanking_dna_seq import selected_frag_dict, add_5prime_end, add_3prime_end
from translate import translate

def get_5_prseq_and_info(row):
    five_side_seqs = []
    frames_and_starts = []
    start_frame = row["qframe"]
    start_codon_pr = row["start_codon_on_sadj_pr"]
    if start_codon_pr != 0:
        return [(start_frame, start_codon_pr)], ['']
    
    start_codon_pr = 1 
    frame = 1
    frames_and_starts.append((frame, start_codon_pr))
    frame_type = "ori"
    if start_frame == 1: frame_type = "new";
    five_side_pr= translate(add_5prime_end(row, frame_type= "ori"))
    five_side_seqs.append(five_side_pr)
    
    if start_frame != 1:
        five_side_pr= translate(add_5prime_end(row, frame_type= "new"))
        five_side_seqs.append(five_side_pr)
        frames_and_starts.append((start_frame, start_codon_pr))
    return frames_and_starts, five_side_seqs


def get_3_prseq_and_info(row):
    three_side_seqs = []
    frames_and_stops = []
    end_frame = row["qframe_end"]
    stop_codon_pr = row["stop_codon_on_sadj_pr"]
    if stop_codon_pr != 0:
        return [(end_frame, stop_codon_pr)], ['']
    
    three_side_pr = translate(add_3prime_end(row, frame_type= "ori"))
    three_side_seqs.append(three_side_pr)
    stop_codon_pr = len(three_side_pr)
    frame = 1
    frames_and_stops.append((frame, stop_codon_pr))
    
    if end_frame != 1:
        three_side_pr= translate(add_3prime_end(row, frame_type= "new"))
        three_side_seqs.append(three_side_pr)
        frames_and_stops.append((end_frame, len(three_side_pr)))
    return frames_and_stops, three_side_seqs