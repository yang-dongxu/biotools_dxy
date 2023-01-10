import pandas as pd
import bioframe as bf 
import pyranges as pr 


def pyranges_to_bioframe(pr:pr.PyRanges):
    return pr.df.rename(columns = {"Chromosome": "chrom", "Start": "start", "End": "end", "Strand": "strand"})

def bioframe_to_pyranges(bf:pd.DataFrame):
    return pr.PyRanges(bf)
