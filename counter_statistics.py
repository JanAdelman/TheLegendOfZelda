import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr
import numpy as np 
import seaborn as sns 
import matplotlib
import matplotlib.pyplot as plt
from numpy import genfromtxt
import seaborn as sns

path = str(pathlib.Path(__file__).parent.absolute())

def get_normalised_df(df):
    # sum of whole df 
    sum = df.to_numpy().astype('uint64').sum()
    print(sum)
    # sum of the rows 
    row_sum = df.sum(axis = 0, skipna = True)
    
    # get normalised df
    assert (row_sum/sum).sum() == 1

    return row_sum/sum 

# load dataframe for enhancer_regions
counts_enhancer_region = pd.read_csv(path+ '/counts_enhancer_regions_dm3.txt', delimiter=',').set_index('identity')
counts_enhancer_region = counts_enhancer_region.iloc[:,0:16]
# get the normalised counts 
normalised_enhancer_regions = get_normalised_df(counts_enhancer_region)
normalised_enhancer_regions = normalised_enhancer_regions.to_frame()

# load data for the non enhancer regions 
counts_non_enhancer_region = pd.read_csv(path+ '/counts_compelement_enhancer_regions_dm3.txt', delimiter=',').set_index('identity')
counts_non_enhancer_region = counts_non_enhancer_region.iloc[:,0:16]
#get thee normalised counts
normalised_complement_regions = get_normalised_df(counts_non_enhancer_region)
normalised_complement_regions = normalised_complement_regions.to_frame()

# load the counts data for the whole genome 
counts = pd.read_csv(path+'/counts.txt', delimiter=',').set_index('identity')
counts = counts.iloc[:,0:16]
# get the normalised counts for the whole genome 
normalised_whole_genome = get_normalised_df(counts)
normalised_whole_genome = normalised_whole_genome.to_frame()

# plot the results 
plt.figure()
sns.set_theme(style="whitegrid")
plt.suptitle('Dinucleotide frequencies D. Melanogaster dm3', size=22)

plt.subplot(1,3,1,)
sns.barplot(x=normalised_enhancer_regions.index, y=normalised_enhancer_regions.iloc[:,0], data=normalised_enhancer_regions)
plt.ylim(0, 0.11)
plt.title('Enhancer regions')
plt.ylabel('frequency')
plt.xlabel('Dinucleotides')

plt.subplot(1,3,2)
sns.barplot(x=normalised_complement_regions.index, y=normalised_complement_regions.iloc[:,0], data=normalised_complement_regions)
plt.ylim(0, 0.11)
plt.title('Non Enhancer regions')
plt.xlabel('Dinucleotides')

plt.subplot(1,3,3)
sns.barplot(x=normalised_whole_genome.index, y=normalised_whole_genome.iloc[:,0], data=normalised_whole_genome)
plt.ylim(0, 0.11)
plt.title('Whole genome')
plt.xlabel('Dinucleotides')

plt.show()

