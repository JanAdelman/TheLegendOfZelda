import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr
import numpy as np 
import seaborn as sns 
import matplotlib
import matplotlib.pyplot as plt

path = str(pathlib.Path(__file__).parent.absolute())

def create_dict(data):
    '''
    this function creates a dictionary for plotting. 
    keys = different developmental stages
    values = counts of annotated elements in the dataset 
    '''

    number_of_values = data.iloc[:,7:].count()
    index = number_of_values.index
    keys_list = list(index)
    values_list = number_of_values.to_numpy()
    zip_iterator = zip(keys_list, values_list)

    return dict(zip_iterator)
#=========================================================================#
# Get data 
#=========================================================================#

# All the zelda motif hits found in the drosophila genome 
zelda_hits_full = pr.read_bed(str(path) + '/TallyNumberOfSites/pwmscan_dm3_5002_7753.bed', as_df=True)

# data provided by the starch lab of known enhancer elements of drosophila, https://enhancers.starklab.org
enhancer_info = pd.read_csv(str(path) + '/tiles.csv')

# dataset that found the enhancer locations and the zelda locations that match up 
enhancer_zelda_matches = pd.read_csv(str(path) + '/enhander_zelda_matches')


# dictionary of annotated counts and developmental stages for the full dataset
dict_full_dataset = create_dict(enhancer_info)

colnames = list(enhancer_info.columns.values) 
full_matched_df = pd.DataFrame(columns = colnames)

# index the full dataset with the start location of the found enhancers. This allows to recover the whole
# row which can be used to visualise the number of enhancer with known function 
for i in range(0, enhancer_zelda_matches.shape[0]):
    full_matched_df = full_matched_df.append(enhancer_info.loc[enhancer_info['Start'] == enhancer_zelda_matches.iloc[i,1]])

full_matched_df.to_csv(path + '/full_matched_dataframe.csv',sep="\t", index=False)

# dictionary of annotated counts and developmental stages for the full dataset
dictionary = create_dict(full_matched_df)

#=========================================================================#
# Plot the distribution of the annotated enhancers 
#=========================================================================#

plt.bar(range(len(dictionary)), list(dictionary.values()), align='center')
plt.xticks(range(len(dictionary)), list(dictionary.keys()))
plt.title('Distribution of annotated enhancers for zelda subset')
plt.xlabel('developmental stages')
plt.ylabel('known action of enhancer')
plt.show()

plt.bar(range(len(dict_full_dataset)), list(dict_full_dataset.values()), align='center')
plt.xticks(range(len(dict_full_dataset)), list(dict_full_dataset.keys()))
plt.title('Distribution of annotated enhancers for complete enhancer dataset')
plt.xlabel('developmental stages')
plt.ylabel('known action of enhancer')
plt.show()