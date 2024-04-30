def overlap(df, range1cols, range2cols):
    maxstart = df[[range1cols[0], range2cols[0]]].max(axis=1)
    minend   = df[[range1cols[1], range2cols[1]]].min(axis=1)
    return minend - maxstart + 1