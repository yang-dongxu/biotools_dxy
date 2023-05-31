import numpy as np
import pandas as pd
import bioframe as bf 
from bioframe import *
import logging
from functools import partial, wraps
import inspect

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


class BioFrameMethods():
    methods_for_overwrite = [
        "is_bedframe",
        "assign_view","closest","cluster","complement","count_overlaps","coverage"
        "expand","merge","overlap","select","select_indices","select_labels","select_mask",
        "setdiff","sort_bedframe","subtract","trim"
        ]
    params_check = ["chrom_col", "start_col", "end_col", "strand_col"]


    def __init__(self,parent) -> None:
        self.parent = parent
        self.info = dict(
            chrom_col = self.parent._chrom_col, 
            start_col = self.parent._start_col, 
            end_col = self.parent._end_col , 
            strand_col = self.parent._strand_col
            )
        self.cols = (self.parent._chrom_col, self.parent._start_col, self.parent._end_col)

        methods_list = [method_name for method_name in dir(bf) if callable(getattr(bf, method_name))]
        for method_name in methods_list:
            if method_name in self.methods_for_overwrite:
                setattr(self, method_name, self.wrap_params(getattr(bf, method_name)))
            else:
                setattr(self, method_name, getattr(bf, method_name))


    def wrap_cols(self,option="cols",**kwargs):
        params_dict = {option:self.cols}
        params_dict.update(kwargs)
        return params_dict

    def wrap_infos(self,**kwargs):
        params_dict = self.info.copy()
        params_dict.update(kwargs)
        return params_dict
    

    def wrap_params(self,func, *args, **kwargs):
        sig = inspect.signature(func)
        params = sig.parameters
        if len( set(self.params_check) & set(params.keys()) ) != 0:
            new_func =  partial(func, self.parent, *args, **self.wrap_infos(**kwargs))
        else:
            if "cols1" in params.keys():
                new_func =  partial(func, self.parent, *args, **self.wrap_cols(option="cols1", **kwargs))
            elif "cols" in params.keys():
                new_func =  partial(func, self.parent, *args, **self.wrap_cols(option="cols",**kwargs))
        new_func.__doc__ = func.__doc__
        return new_func
            
         
    def resize(self, size, fix = "center", *args, **kwargs):
        return self.wrap_params(resize)(size, fix, *args, **kwargs)

    


class BioRegions(pd.DataFrame):
    params_kept = ["_chrom_col", "_start_col", "_end_col", "_strand_col"]
    params_kept_default = ["chrom", "start", "end", "strand"]

    def __init__(self,*args,**kwargs) -> None:
        raw_params = kwargs.copy()

        for param in BioRegions.params_kept:
            if param in kwargs.keys():
                del kwargs[param]
        super().__init__(*args,**kwargs)

        for param, pdault in zip(BioRegions.params_kept, BioRegions.params_kept_default):
            self.__setattr__(param, pdault)
            if param in raw_params.keys():
                setattr(self, param, kwargs[param])


        for param in BioRegions.params_kept:
            assert getattr(self, param) in self.columns, f"{param=} not in columns! The column name should be set by {param} = 'column_name'"

        self.bf = BioFrameMethods(self)

@wraps
def wrap_bf(func):
    def wrapper(*args, **kwargs):
        return BioRegions(func(*args, **kwargs))
    return wrapper

read_table = wrap_bf(bf.read_table)
