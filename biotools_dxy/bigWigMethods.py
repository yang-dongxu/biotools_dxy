
import numpy as np
import pandas as pd

import bbi 
# import pyBigWig as pbw

from matplotlib import pyplot as plt


def __bw_getSignal_bins_kernal(bw, regions, chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand", strand = True, scale = False, engine = "bbi", bins = 1, missing = 0):
    if engine not in ["pybbi","bbi", "pbw", "pybigwig"]:
        raise ValueError(
            "engine should be one of ['pybbi','bbi', 'pbw', 'pybigwig'], but got %s" % engine
            )
    if engine in ["pybbi","bbi"]:
        with bbi.open(str(bw)) as bwf:
            mtx = bwf.stackup(regions[chrom_col],regions[start_col],regions[end_col], bins=bins, missing=missing)
            if scale == True:
                mean_ = bwf.info["summary"]["mean"]
                mtx = mtx/mean_
    elif engine in ["pbw", "pybigwig"]:
        raise NotImplementedError("pybigwig is not implemented yet")
        with pbw.open(str(bw)) as bwf:
            values = []
            for c, s, e in zip(regions[chrom_col],regions[start_col],regions[end_col]):
                values.append(bwf.stats(c,s,e, nBins = bins, missing = missing))
            mtx = np.array(values)
            if scale == True:
                mean_ = bwf.header()["sumData"]/bwf.header()["nBasesCovered"]
                mtx = mtx/mean_
    if strand_col in regions.columns and strand == True:
        print("using strand...")
        mtx[regions[strand_col] == "-",:] = mtx[regions[strand_col] == "-",::-1]
    
    return mtx


def bw_getSignal_bins(
    bw, regions:pd.DataFrame, 
    chrom_col = "chrom",start_col = "start", end_col = "end", strand_col = "strand", 
    width = 4000,bins = 100, missing = 0, 
    scale = False, strand = True, 
    inherit_cols = ["name"],
    prefix = "bin_", engine = "bbi"
    ):
    assert chrom_col in regions.columns
    assert start_col in regions.columns
    assert end_col in regions.columns
    regions_raw = regions
    regions = regions.copy()
    if width > 0:
        mid = (regions[start_col] + regions[end_col])//2
        regions[start_col] = mid - int(width/2)
        regions[end_col] = mid + int(width/2)

        mtx = __bw_getSignal_bins_kernal(bw, regions, 
            chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col, 
            strand = strand, scale = scale, bins = bins, missing = missing, engine = engine
            )
        
        allCols = [f"{prefix}{i}" for i in range(mtx.shape[1])]
        df_signal = pd.DataFrame(data = mtx, columns = allCols)

        # if inherit_cols is "all": then inherit all columns
        if inherit_cols == "all":
            inherit_cols = regions_raw.columns
        used_cols = [x for x in inherit_cols if x in regions_raw.columns]
        for c in used_cols:
            try:
                df_signal[c] = regions_raw[c].values
            except Exception as e:
                print(f"Error: {c} not in regions")
                raise(e)
    return df_signal[used_cols + allCols]


def plotHeatmapSplit(df,signalCols,splitCol="label",splitOrder =None, splitLabels = None,orderCol="order", orderAscending = False, fig = None, xlabel = "Distance from motif center (bp)", xticks = None, xticklabels = ["-2k","0","2k"], default_figsize = (3,9), *args, **kwargs):
    if fig is None:
        fig = plt.figure(figsize = default_figsize)
    dfg = df.groupby(splitCol)
    sizes = dfg.size().to_dict()
    if splitOrder is None:
        splitOrder = list(sizes.keys() )
        splitOrder = [i for i in splitOrder if sizes[i] > 0] # drop blank
    if splitLabels is None:
        splitLabels = splitOrder
    else:
        assert len(splitOrder) == len(splitLabels)
        print("use self-defined labels ...")
    splitLabelsDict = {i:j for i,j in zip(splitOrder, splitLabels)}
    splitOrder = [i for i in splitOrder if i in sizes.keys()]
    heightRatios = [sizes[x] for x in splitOrder]

    axs = fig.subplots(len(heightRatios), 1, gridspec_kw={'height_ratios': heightRatios})
    for i,ax in enumerate(axs):
        cluster = splitOrder[i]
        dfc = dfg.get_group(cluster)

        if isinstance(orderCol, str):
            if orderCol in dfc.columns:
                dfc = dfc.sort_values(orderCol,ascending = orderAscending)
            else:
                print(f"warning: {orderCol} not in dfc.columns, skip sorting")
        elif isinstance(orderCol, list):
            assert all([x in dfc.columns for x in orderCol])
            dfc = dfc.sort_values(orderCol,ascending = orderAscending)
        else:
            print(f"error: oderCol should be str or list, but got {type(orderCol)}")
            raise ValueError
        ax.imshow(dfc[signalCols], aspect="auto", *args, **kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(splitLabelsDict[cluster])
    # if xticks is none, turn it to [0,mid, end]
    if xticks is None:
        xticks = [0, int(len(signalCols)/2), len(signalCols)-1]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel(xlabel)
    return fig, ax    



