#! /usr/bin/env python 

'''
Anneliek ter Horst
python script to get only the <75% coverage contigs from the bed file csv
bedtools genomecov outputs a tsv that is the input file here. This script 
obtains the average coverage per contig. If its more than 75%, we will discard the contig.
'''

# imports
from __future__ import division
import sys
import pandas as pd

# open the file to read with pandas
df = pd.read_csv('infile.tsv',sep='\t', header=None)

# add column names so we know what we are looking at
df.columns = ['sequence', 'start_pos', 'stop_pos' , 'coverage']

# Add a column with the number of nucleotides that share the same coverage
df['num_of_nucl'] = df['stop_pos'] - df['start_pos']

# for each sequence, add a extra column witht the sequences total length. Calculating percentages
# becomes easier this way
df['total_len_seq'] = df.groupby(['sequence'],sort=False)['stop_pos'].transform(max)

# calculate percentage coverage for each mapping piece to the contig
df['percentage_of_seq'] = df['num_of_nucl'] / df['total_len_seq']


# make a new df where only values with coverage not 0 are in
df_no0 = df.loc[df['coverage'] != 0]

# add up all the percentages for sequences with the same name
# reset index to keep column names
df_percentages = df_no0.groupby(['sequence'])['percentage_of_seq'].agg('sum').reset_index()

# make it actual percentages by multiplying with 100
df_percentages['percentage_of_seq'] = df_percentages['percentage_of_seq']*100


# Select the cases where the percentage coverage is < 74.99
percentage_under_75 = df_percentages['percentage_of_seq'] < 74.99

# keep cases where perc coverage < 75, to new csv
df_percentages[percentage_under_75].to_csv(sys.argv[2], sep='\t')
