#!/bin/python

import pandas as pd
import glob
import os
#import functools

if not os.path.exists("data"):
    os.mkdir("data")

if os.path.exists("snakemake/ASSEMBLY/ISOQUANT"):
#glob all count files
    allFiles = glob.glob("snakemake/ASSEMBLY/ISOQUANT/*/*.gene_counts.tsv")
    combined = pd.concat([pd.read_csv(f, sep="\t", header=1, usecols=[0,1], names=['feature_id', str(f).replace('.txt', '').split('/')[4]]) for f in allFiles], axis=1)
    combined = combined.loc[:,~combined.columns.duplicated()] #remove duplicate columns
    combined.to_csv("data/isoquant_gene_counts.tsv", index=False, header=True, sep='\t') #output