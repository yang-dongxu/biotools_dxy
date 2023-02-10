import numpy as np
import pandas as pd
import bioframe as bf 
import logging

def __get_fix_position(
    df: pd.DataFrame,
    fix = "center",
    chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand",
    ):
    '''
    get the fixed position for the genomic regions
    df: a dataframe with genomic regions
    fix: one of "center", "left", "right", "start", "end".
        "start" and "end" means the start and end of the original region according to the strand, so strand_col is required
    '''
    for col in [chrom_col, start_col, end_col]:
        if col not in df.columns:
            raise ValueError(f"{col} not in columns")
    if fix in ["start", "end"] and strand_col not in df.columns:
        raise ValueError(f"{strand_col} not in columns")
    df = df.copy()
    position = None
    if fix == "center":
        position = (df[start_col] + df[end_col])//2
    elif fix == "left":
        position = df[start_col]
    elif fix == "right":
        position = df[end_col]
    elif fix == "start":
        position = df.apply(lambda x: x[start_col] if x[strand_col] == "+" else x[end_col], axis = 1)
    elif fix == "end":
        position = df.apply(lambda x: x[end_col] if x[strand_col] == "+" else x[start_col], axis = 1)
    else:
        raise ValueError(f"fix {fix} not supported")
    return position.to_numpy()


def __generate_new_regions(position:np.array, size, strands = None, fix = "center"):
    if fix == "center":
        start = position - size//2
        end = position + size//2
    elif fix == "left":
        start = position
        end = position + size
    elif fix == "right":
        start = position - size
        end = position
    elif fix in [ "start", "end"]:
        if len(strands) != len(position):
            raise ValueError("uneuqal length of input position and strand!")
        unique_strand = set(strands)
        if len(unique_strand.difference(set(["+","-"]))) != 0:
            logging.getLogger(__name__).warning(f"strands contains {unique_strand}!")
        v1 = position
        shifts = np.repeat(size, v1.size)
        strand = np.array([-1 if s == "-" else 1 for s in strands])
        shifts = shifts*strand
        if fix == "start":
            v2 = v1 + shifts
        else:
            v2 = v1 - shifts
        start = np.minimum(v1, v2)
        end = np.maximum(v1,v2)
    else:
        raise ValueError(f"fix {fix} not supported")
    return start, end





# method like resize of GenomicRanges in R
def resize(
    df: pd.DataFrame, 
    size: int, fix = "center", 
    chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand",
    *args, **kwargs
    ):
    '''
    resize the genomic regions to a new size
    df: a dataframe with genomic regions
    size: the new size for the regions, result may be different according to the method
    fix: one of "center", "left", "right", "start", "end". 
        "start" and "end" means the start and end of the original region according to the strand, so strand_col is required
    '''
    for col in [chrom_col, start_col, end_col]:
        if col not in df.columns:
            raise ValueError(f"{col} not in columns")
    if fix in ["start", "end"] and strand_col not in df.columns:
        raise ValueError(f"{strand_col} not in columns")
    df = df.copy()
    position = __get_fix_position(df, fix, chrom_col, start_col, end_col, strand_col)
    if strand_col in df.columns:
        strands = df[strand_col]
    else:
        strands = None
    start, end = __generate_new_regions(position, size = size, strands = strands, fix = fix)
    df[start_col] = start
    df[end_col] = end

    return df

