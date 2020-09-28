#! /usr/bin/env python 
# imports
import pandas as pd
import numpy as np
import re
import sys


# open the coverage table
df = pd.read_csv('infile.tsv' ,sep='\t')

# only keep entries with a virfinder score over 0.9
df = df[df.score > 0.899]
# only keep entries with a p value under 0.05
df = df[df.pvalue < 0.05]

# print how many entries are left after those tresholds
print len(df)

# save results in new file
df.to_csv('outfile.tsv', sep='\t')