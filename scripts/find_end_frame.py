def find_end_frame(qseq, start_frame):
    fs_num = qseq.count("/")
    bs_num = qseq.count("\\")
    end_frame = (start_frame - fs_num + bs_num)%3
    if end_frame==0:
        return 3
    return end_frame