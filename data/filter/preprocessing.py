########################################################

# Code and data for 
# ‘Temperature-driven scaling patterns emerge in diatom 
# gene expression across the global ocean'

########################################################

# Python script that transforms the original meta-omics 
# abundance data into abundance data per sample. 

########################################################

# Required modules:
import pandas as pd
import dask.dataframe as dd
import os

########################################################

# Create the directories for the new abundance tables
# per sample.
os.makedirs('metaG_micro/by_sampleID', exist_ok=True)
os.makedirs('metaT_micro/by_sampleID', exist_ok=True)

# Read the original abundance data.
# list of files in the ssv directory that are .ssv files
# the original data was cut into pieces to make them more
# manageable
ssv_files = [f for f in os.listdir('ssv') if f.endswith('.ssv')]
# sort the list of files
ssv_files.sort()

for file in ssv_files[:]:
    # read the subfile into a dataframe
    df = dd.read_csv('ssv/' + file, sep=' ', header=0)
    # get the unique sample IDs
    unique_values = df.iloc[:,1].unique().compute()
    
    # for each sample in the subfile
    for index in range(len(unique_values)):

        value = unique_values[index]
        # extract all genes in this sample
        df2 = df.loc[df.iloc[:, 1] == value, [df.columns[0], df.columns[2]]]
        # convert dask dataframe to pandas dataframe
        df2 = df2.compute()

        # correctly assign the sample to metaG and metaT 
        if ((value[-2:] == '11') or (value[-2:] == '12')):   
            filename = 'metaG_micro/by_sampleID/' + value + '.csv'
        elif ((value[-2:] == '13') or (value[-2:] == '14')):
            filename = 'metaT_micro/by_sampleID/' + value + '.csv'
        # if the file already exists, append to it. 
        # otherwise, create a new file
        # because some samples are divided over two ssv subfiles
        if os.path.isfile(filename):
            df2.to_csv(filename, sep=',', header=False, mode='a', index=False)
        else:
            df2.to_csv(filename, sep=',', header=True, index=False)