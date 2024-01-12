import os 
import sys 
import pandas as pd 
import pyranges as pr
import bioframe as bf
from .utilities import regions

import logging

def parse_refGTF(ref_gtf, chrom_size, tss_upstream = 2000, tss_downstream = 2000, tes_upstream = 2000, tes_downstream = 2000, distal = 10_000, *args, **kwargs):
    """
    Parse reference annotation file in gtf format.
    Parameters
    ----------
    ref_gtf : str
        Reference annotation file in gtf/gtf.gz format. 
    tss_upstream : int
        Upstream region of TSS to annotate.
    tss_downstream : int
        Downstream region of TSS to annotate.
    Returns
    -------
    bioframe: promoter, exon, intron, intergenic, distal, and pyranges:ref_gtf in a dict format
        
    """
    if not os.path.exists(ref_gtf):
        raise(FileNotFoundError("Reference annotation file does not exist!"))
    
    logging.getLogger(__name__).info("Parsing reference annotation file...")
    pr_refgtf = pr.read_gtf(ref_gtf)
    assert isinstance(pr_refgtf, pr.PyRanges), "Reference annotation file is not in gtf format!"
    # check chromsome format. if not chr* format, add chr
    chromosomes = pr_refgtf.chromosomes
    if not all([i.startswith("chr") for i in  chromosomes]):
        logging.getLogger(__name__).warning("Formating reference gtf chromosome names...")
        df = pr_refgtf.df
        df["Chromosome"] = "chr" + df["Chromosome"].astype(str)
        pr_refgtf = pr.PyRanges(df)


    # parse chromszie 
    if chrom_size in ["mm9","mm10","hg19","hg38"]:
        chromSize = bf.fetch_chromsizes(chrom_size)
    else:
        assert os.path.exists(chrom_size), "chromSize file does not exist!"
        chromSize = bf.read_chromsizes(chrom_size)
    # to pyranges
    pr_chromSize = pr.PyRanges(
        chromSize.to_frame().assign(Chromosome = chromSize.index, Start = 0, End = chromSize.values).reset_index(drop = True)
        )

    # get features: tss, exon, intron, intergenic
    pr_tss = pr_refgtf.features.tss()
    pr_tes = pr_refgtf.features.tes()
    pr_exon = pr_refgtf[pr_refgtf.Feature == "exon"]
    pr_refgtf.features.introns(by = "transcript")
    pr_intergenic = pr_chromSize.subtract(pr_refgtf)

    # get promoter from tss and distal from intergenic
    pr_promoter = pr_tss.extend({"5": tss_upstream, "3": tss_downstream})
    pr_tes = pr_tes.extend({"5": tes_upstream, "3": tes_downstream})
    pr_distal = pr_intergenic[pr_intergenic.lengths() > 2*distal].extend(-distal) # shrink to distal from gene

    # turn pyranges to bioframe 
    bf_exon = regions.pyranges_to_bioframe(pr_exon)
    bf_intron = regions.pyranges_to_bioframe(pr_intron)
    bf_intergenic = regions.pyranges_to_bioframe(pr_intergenic)
    bf_promoter = regions.pyranges_to_bioframe(pr_promoter)
    bf_tes = regions.pyranges_to_bioframe(pr_tes)
    bf_distal = regions.pyranges_to_bioframe(pr_distal)

    out = {
        "promoter": bf_promoter,
        "tes": bf_tes,
        "exon": bf_exon,
        "intron": bf_intron,
        "intergenic": bf_intergenic,
        "distal": bf_distal,
        "ref": pr_refgtf,
    }
    return out


def byGene(
    feature, ref_gtf, chrom_size = "mm10", tss_upstream = 2000, tss_downstream = 2000, tes_upstream = 2000, tes_downstream = 2000, distal = 10_000, 
    chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand",
    use_strand:bool = False, count:bool = False, *args, **kwargs
):
    """
    Annotate genomic distribution of gene-related features. 
    Parameters
    ----------
    feature : pandas.DataFrame
        A dataframe with at least a column named 'chrom' and 'start' and 'end', or given by chrom_col, start_col, end_col.
    ref_gtf : str
        Reference annotation file in gtf/gtf.gz format.
    tss_upstream : int
        Upstream region of TSS to annotate.
    tss_downstream : int
        Downstream region of TSS to annotate.
    tes_upstream : int
        Upstream region of TES to annotate.
    tes_downstream : int
        Downstream region of TES to annotate.
    chrom_col : str
        Column name of chromosome.
    start_col : str
        Column name of start position.
    end_col : str
        Column name of end position.
    strand_col : str
        Column name of strand.
    use_strand : bool
        Consider strand when annotating.
    count : bool
        Whether to count the number of features in each region, or just True/False.
    Returns
    -------
    feature : pandas.DataFrame
        Input DataFrame with additional columns: 'promoter', 'exon', 'intron', 'intergenic', 'distal'.
    """
    # parse ref_gtf
    refs = parse_refGTF(ref_gtf, chrom_size = chrom_size, tss_upstream = tss_upstream, tss_downstream = tss_downstream, tes_upstream = tes_upstream, tes_downstream = tes_downstream, distal = distal, *args, **kwargs)

    feature = feature.copy()

    onterm = None
    if use_strand:
        onterm = [strand_col]
    # get overlaps by bf.count_overlaps
    for region in ["promoter","tes", "exon", "intron", "intergenic", "distal"]:
        feature[region] = bf.count_overlaps(
            feature, refs[region],
            cols1=(chrom_col, start_col, end_col), on = onterm, return_input=False
        )["count"]
        if not count:
            feature[region] = feature[region].astype(bool)
    return feature

def stat_features(regions: pd.DataFrame, features = ["promoter","tes", "exon", "intron", "distal"], by=[]):
    selected_cols = [col for col in features if col in regions.columns]
    grouped_cols = [col for col in by if col in regions.columns]

    if len(selected_cols) == 0 :
        raise(ValueError("No columns selected in the given table!"))
    
    def stat(regions, features):
        out = {}
        all_num = len(regions)
        region_selected = regions.copy()
        for col in features:
            region_selected = region_selected[region_selected[col] > 0 ]
            out[col] = all_num - len(region_selected)
            all_num = len(region_selected)
        out["other"] = all_num
        return out
    
    if len(grouped_cols) == 0:
        dfo = stat(regions, features)
        dfo_long = pd.DataFrame(dfo, index = [0]).stack().to_frame().reset_index().rename(columns = {"level_1": "feature", 0: "count"}).drop(columns="level_0")
        dfo_long["percent"] = dfo_long["count"]/dfo_long["count"].sum()
        
    else:
        outs = []
        for g,df in regions.groupby(grouped_cols):
            dict_out = stat(df, features)
            if isinstance(g, tuple):
                for col, val in zip(grouped_cols, g):
                    dict_out[col] = val
            elif isinstance(g, str):
                dict_out[grouped_cols[0]] = g
            else:
                raise(TypeError("Unknown grouped columns type!"))
            outs.append(dict_out)
        dfo = pd.DataFrame(outs)
        dfo_long = dfo.set_index(grouped_cols).stack().to_frame()\
            .reset_index().rename(columns = {"level_2": "feature", 0: "count"})
        dfo_long["percent"] = dfo_long.groupby(grouped_cols,group_keys=False)["count"].apply(lambda x: x/x.sum())
    dfo_long["feature"] = pd.Categorical(dfo_long["feature"], features + ["other"])
    return dfo, dfo_long


    


def parse_repeats(
    repeats_txt, 
    repeats_class = ["LINE","SINE","LTR"], 
    repeats_families = ["Alu", "L1"], 
    repeats_subfamilies = ["B3A"], *args, **kwargs
):
    """
    Parse repeats annotation file in txt format.
    Parameters
    ----------
    repeats_txt : str
        Repeats annotation file in txt format, can be downloaded from UCSC table browser (rmsk).
    repeats_class : list
        Repeats class to annotate
    repeats_families : list
        Repeats families to annotate
    repeats_subfamilies : list
        Repeats subfamilies to annotate
    Returns
    -------
    A dictionary of repeats annotation. for class, family, subfamily, key is c/f/s_termname, value is a bioframe with strand
    """
    if not os.path.exists(repeats_txt):
        raise(FileNotFoundError("Repeats annotation file does not exist!"))
    df_repeats = pd.read_csv(repeats_txt, sep = "\t", header = None,usecols=[5,6,7,9,10,11,12])
    df_repeats.columns = ["chrom","start","end","strand","subfamily","class","family"]

    out = {}
    out["a_all"] = df_repeats
    for c in repeats_class:
        out[f"c_{c}"] = df_repeats[df_repeats["class"] == c]
    for f in repeats_families:
        out[f"f_{f}"] = df_repeats[df_repeats["family"] == f]
    for s in repeats_subfamilies:
        out[f"s_{s}"] = df_repeats[df_repeats["subfamily"] == s]
    return out


    

def byRepeats(
    feature, repeats_txt, 
    repeats_class = ["LINE","SINE","LTR"], 
    repeats_families = ["Alu", "L1"], 
    repeats_subfamilies = ["B3A"],
    chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand",
    use_strand:bool = False, count:bool = False, *args, **kwargs

):
    """
    Annotate genomic distribution of features.
    Parameters
    ----------
    feature : pandas.DataFrame
        A dataframe with at least a column named 'chrom' and 'start' and 'end'.
    repeats_txt : str
        Repeats annotation file in txt format.
    repeats_class : list
        Repeats class to annotate
    repeats_families : list
        Repeats families to annotate
    repeats_subfamilies : list
        Repeats subfamilies to annotate
    chrom_col : str
        Column name of chromosome.
    start_col : str
        Column name of start position.
    end_col : str
        Column name of end position.
    strand_col : str
        Column name of strand.
    use_strand : bool
        Whether to use strand information.
    count: bool
        Whether to count the number of overlaps, or just True/False
    Returns
    -------
    feature : pandas.DataFrame
        A dataframe with additional selected repeats columns.
    """
    feature = feature.copy()
    repeats = parse_repeats(repeats_txt, repeats_class, repeats_families, repeats_subfamilies)
    onterm = None
    if use_strand:
        onterm = [strand_col]
    
    for r in repeats:
        feature[r] = bf.count_overlaps(
            feature, repeats[r],
            cols1=(chrom_col, start_col, end_col), on = onterm, return_input=False
        )["count"]
        if not count:
            feature[r] = feature[r].astype(bool)
    return feature




def genomicDistributionAnnotation(
    feature, ref_gtf, repeats_txt, tss_upstream = 2000, tss_downstream = 2000,
    chrom_col = "chrom", start_col = "start", end_col = "end", strand_col = "strand",
    use_strand:bool = False, count:bool = False, *args, **kwargs
):
    """
    Annotate genomic distribution of features.
    Parameters
    ----------
    feature : pandas.DataFrame
        A dataframe with at least a column named 'chrom' and 'start' and 'end'.
    ref_gtf : str
        Reference annotation file in gtf/gtf.gz format.
    repeats_txt : str
        Repeats annotation file in txt format.
    tss_upstream : int
        Upstream region of TSS to annotate.
    tss_downstream : int
        Downstream region of TSS to annotate.
    chrom_col : str
        Column name of chromosome.
    start_col : str
        Column name of start position.
    end_col : str
        Column name of end position.
    strand_col : str
        Column name of strand.
    use_strand : bool
        Whether to use strand information.
    count: bool
        Whether to count the number of overlaps, or just True/False
    Returns
    -------
    featuresGene: pandas.DataFrame
        A dataframe with additional gene columns.
    featuresRepeats: pandas.DataFrame
        A dataframe with additional repeats columns.

    """
    
    # validata params

    if not isinstance(feature, pd.DataFrame):
        raise(TypeError("feature must be a pandas.DataFrame!"))

    for c in [chrom_col, start_col, end_col]:
        if c not in feature.columns:
            raise(ValueError("feature must have a column named '{}'!".format(c)))
    if use_strand == True and strand_col not in feature.columns:
        raise(ValueError("feature must have a column named '{}'!".format(strand_col)))

    features = feature.copy()
    featuresGene = byGene(features, ref_gtf, chrom_col, start_col, end_col, strand_col, use_strand, count)
    featuresRepeats = byRepeats(features, repeats_txt, chrom_col, start_col, end_col, strand_col, use_strand, count)

    return featuresGene, featuresRepeats
