#! /usr/bin/env python 

'''
python script to remove vOTUs from coverage table that do not have >75% 1x coverage.
Infiles are the is created by percentage_coverage.py, the list of vOTUs covered under 75% and 
the coverage table itself. 
Outfile is a coverage table with vOTUs that are all at least 1x covered over 75% of their lenght.
'''

# imports
import pandas as pd
import numpy as np
import re
import sys

# open the coverage table
df = pd.read_csv('coverage_table.tsv',sep='\t')

# open the under 75% coverage list per sample
df_under_75= pd.read_csv('under_75_per_sample.csv',sep='|')

# make the headers readable, cut from the _S# part
# this might need to be adapted or left out depending on your files
df.columns = df.columns.str.split('_L006').str[0] + ''
df_under_75.columns = df_under_75.columns.str.split('_L006').str[0] + ''


#remove viral contigs with all 0s from the coverage file
df_cov = remove_all_zeros(df)


# check the lenght of the dataframe and the number of non zero values in each 
# column (for each sample)
print 'before', len(df_cov), df_cov.astype(bool).sum(axis=0)


#Do the function that checks each column(sample) for a coverage under 75%. 
for col in df_under_75.columns:
       df_cov[col] = df_cov[['#contig', col]].apply(lambda x: check_similarity(x['#contig'], df_under_75[col].values[1:], x[col]),axis=1)

    
#after removing these under 75% coverage ones, remove the columns with all 0 again.
df_cov = remove_all_zeros(df_cov)


# check the len of df and num of viral contigs with > 75% coverage for each sample
print 'after', len(df_cov),  df_cov.astype(bool).sum(axis=0)

#print df_cov
df_cov.to_csv('biochar_covtab_filter.csv', sep='|')


# function to remove all 0s from dataframe. If a vOTU only has zero's, remove it ( so not covered)
def remove_all_zeros(df):

    # get a list of columns that could contain zero
    all_columns =  list(df.columns)[2:]

    # check if all columns contain zero
    df["all_zeros"] = (df[all_columns] == 0).all(axis=1)

    # remove all that have only zeros
    df = df.loc[df["all_zeros"] == False]

    # remove all zeros column
    df.drop(['all_zeros'], axis=1, inplace=True)
    
    return df

# check if the vOTU is in the list of under 75% covered, if so, remove that vOTU from table
def check_similarity(contact_name, contact_names_to_check, curent_value):
        
    if contact_name in contact_names_to_check:
        return 0
    else: return curent_value