# import necessary packages
import numpy as np
import os
import csv
import re 
import glob
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# get working directory
cwd = os.getcwd()
my_path = os.path.abspath(__file__)

def get_TF_data(TF_data):
    '''
    Function to return a pandas data frame for one TF, the data has not the same column number 
    for each row thus column are appended such that the dimensions match 
    '''

    with open(fname, newline='') as temp_f:
        # get # of columns in each line
        col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]
           
        
    column_names = [i for i in range(0, max(col_count))]
    df = pd.read_csv(fname, header=None, delimiter="\t", names=column_names)
    name = str(df[1][1])
    df = df.drop([0,1,2,3,4,5,6])

    df = df.to_numpy()
    df = np.delete(df, 0, 1)
    df = df.astype(np.float)

    return df, name

path = "/Users/janadelmann/polybox/LabRotation1/TF_PWM/*.csv"
# get data for each TF
TF_info = pd.DataFrame(columns=['TF_Name', "Info"])
TF_info_normalised = pd.DataFrame(columns=['TF_Name', 'Info_normalised'])

def update_pvm(pvm_matrix, L, n_nucleotides):
    
    '''
    Ineficient Function to add pseudo probabilities to nucleotides which occur with probability 0 
    in a row. 
    '''

    for l in range(0, L):
        count = np.count_nonzero(pvm_matrix[l][:] != 0)

        if count == 4:
            continue

        if count == 3:
            for nucleotide in range(0, n_nucleotides): 
                if pvm_matrix[l][nucleotide] == 0:
                    pvm_matrix[l][nucleotide] += 10**(-6)
                else:
                    pvm_matrix[l][nucleotide] -= (10**(-6))/3

        if count == 2:
            for nucleotide in range(0, n_nucleotides): 
                if pvm_matrix[l][nucleotide] == 0:
                    pvm_matrix[l][nucleotide] += 10**(-6)
                else:
                    pvm_matrix[l][nucleotide] -= (10**(-6))

        if count == 1:
            for nucleotide in range(0, n_nucleotides): 
                if pvm_matrix[l][nucleotide] == 0:
                    pvm_matrix[l][nucleotide] += 10**(-6)
                else:
                    pvm_matrix[l][nucleotide] -= (3*10**(-6))
                if count == 0:
                    print('alarm')

    return pvm_matrix

# https://www.ncbi.nlm.nih.gov/genome?term=vih&cmd=DetailsSearch
# median GC-content of ID-47, 42%
nucleo_distr = [0.29, 0.21, 0.21, 0.29]

nucleo_unif = 0.25

#==============================================================================#
# Calculate the Information content for each TF
#==============================================================================#

num = 0 

for fname in glob.glob(path):

    filesize = os.path.getsize(fname)
    
    # if file is empty, skip 
    if filesize == 0:
        pass
    else: 
        # get the PVM Matrix + the name of the TF 
        pvm_matrix , name = get_TF_data(fname)

        # print the PVM of Zelda (2 experiments)
        #if name == 'vfl':
            #print(name)
            #print(pvm_matrix)

        L = pvm_matrix.shape[0]
        n_nucleotides = pvm_matrix.shape[1]

        # wacky version on how to add pseudo probabilites if the occurence of a certain nucleotide is 0
        # Remark, log not defined for 0. 
    
        if num in pvm_matrix:
            pvm_matrix = update_pvm(pvm_matrix, L, n_nucleotides)
      
        # calculate the information content of each TF 

        info_tf = 0
        for l in range(0, L):
            sub_sum = 0
            for nucleotide in range(0, n_nucleotides):
                if pvm_matrix[l][nucleotide] == 0:
                    pass 
                else:
                    sub_sum += pvm_matrix[l][nucleotide]*np.log2(pvm_matrix[l][nucleotide]/nucleo_distr[nucleotide])

            info_tf += sub_sum

        
        df = pd.DataFrame([[name, info_tf]], columns = ['TF_Name', "Info"])
        TF_info = TF_info.append(df, ignore_index=True)
        TF_info = TF_info.sort_values(by='TF_Name')
        TF_info.to_csv('/Users/janadelmann/polybox/LabRotation1/information_content_TF', index=False)

        # create dataframe normailsed by length of each motfi 
        info_tf_normalised = info_tf/L
        df_normalised = pd.DataFrame([[name, info_tf_normalised]], columns = ['TF_Name', 'Info_normalised'])
        TF_info_normalised = TF_info_normalised.append(df_normalised, ignore_index=True)
        TF_info_normalised = TF_info_normalised.sort_values(by='TF_Name')
        TF_info_normalised.to_csv('/Users/janadelmann/polybox/LabRotation1/information_content_TF_normalised', index=False)

#==============================================================================#
# Calculate basic statistics from the data TF
#==============================================================================#

# get the mean over all epxeriments 
mean_data = TF_info.groupby(by="TF_Name").mean()
mean_data_normalised = TF_info_normalised.groupby(by='TF_Name').mean()
# get std of mean over experiments 
std_data = mean_data.std()

# finding TFs with big and small info conten 
subset_normalised = mean_data_normalised.sort_values(by='Info_normalised')
subset = mean_data.sort_values(by='Info')
print('=====================')
print('Min/Max info content for normalised TFs')
print(subset_normalised[0:5])
print(subset_normalised[-5:])
print('normalised vfl: ')
print(mean_data_normalised['Info_normalised']['vfl'])
print('====================')
print('Min/Max info content for TFs')
print(subset[0:5])
print(subset[-5:])
print('vfl: ')
print(mean_data['Info']['vfl'])

# test for normality of the data (data is not normally distributed)
shapiro_test_full = stats.shapiro(TF_info['Info'])
shapiro_test_mean = stats.shapiro(mean_data)

# get the experiment with min informataion
min_data = TF_info.groupby(by="TF_Name").min()

# get the experiment with max information 
max_data = TF_info.groupby(by="TF_Name").max()

# Calculate average over the averaged dataset 
mean_of_mean_dataset = mean_data['Info'].mean()
print(mean_of_mean_dataset)
mean_of_normalised_datatset = mean_data_normalised['Info_normalised'].mean()
print(mean_of_normalised_datatset)

# Calulate average over full dataset with repeated experiments for some TFs
mean = TF_info['Info'].mean()
print(mean)

# get information value for Zelda transcription factor 
# mean would be the best 

mean_vfl = mean_data['Info']['vfl']
max_vfl = max_data['Info']['vfl']
min_vfl = min_data['Info']['vfl']
normalised_vfl = mean_data_normalised['Info_normalised']['vfl']

# get information value for GrainyHead transcription factor 
mean_grh = mean_data['Info']['grh']
max_grh = max_data['Info']['grh']
min_grh = min_data['Info']['grh']
normalised_grh = mean_data_normalised['Info_normalised']['grh']

# what fraction of TF are below zelda? 
def sum(df):
    return (df.Info < df['Info']['vfl']).sum()

def sum_normalised(df):
    return (df.Info_normalised < df['Info_normalised']['vfl']).sum()

smaller_zelda = sum(mean_data)
smaller_zelda_normalised = sum_normalised(mean_data_normalised)

# print the size of the DF and print the number of TFs with smaller info content than zelda 
print('=======================================')
print('TFs smaller infoconnntent than Zelda: ')
print(mean_data.shape[0])
print(smaller_zelda)

print('========================================')
print('Normalised TFs smaller infoconnntent than Zelda: ')
print(mean_data_normalised.shape[0])
print(smaller_zelda_normalised)

#==============================================================================#
# Histogram Plot for different dataset compositions
#==============================================================================#


# Plot of only average information 
sns.histplot(mean_data['Info'], alpha = 0.5)
plt.axvline(mean_of_mean_dataset, color='red')
plt.axvline(mean_vfl, color='darkblue')
plt.axvline(max_grh, color='green', alpha = 1)
plt.annotate('Zelda Information Content', xy=(float(mean_vfl), 20), xytext=(float(mean_vfl)+2, 40), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Grainy Head Information Content', xy=(float(mean_grh), 25), xytext=(float(mean_grh)+4.3, 55), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Mean', xy=(float(mean_of_mean_dataset), 15), xytext=(float(mean_grh)+8, 20), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.ylabel("Count",fontsize=14)
plt.xlabel('Information content (bits)', fontsize = 14)
plt.title('Mean information content of D. Melanogaster TFs', fontsize = 17)
plt.savefig('/Users/janadelmann/polybox/LabRotation1/figures/mean_info_content.png')
plt.show()

# plot of the normalised information content 
sns.histplot(mean_data_normalised['Info_normalised'], alpha = 0.5)
plt.axvline(mean_of_normalised_datatset, color='red')
plt.axvline(normalised_vfl, color='darkblue')
plt.axvline(normalised_grh, color='green', alpha = 1)
plt.annotate('Zelda Information Content', xy=(float(normalised_vfl), 20), xytext=(float(normalised_vfl)+0.1, 45), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Grainy Head', xy=(float(normalised_grh), 25), xytext=(float(normalised_grh)-0.5, 55), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Mean', xy=(float(mean_of_normalised_datatset), 15), xytext=(float(normalised_grh)-0.5, 20), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.ylabel("Count",fontsize=14)
plt.xlabel('Information content (bits)', fontsize = 14)
plt.title('Mean information content of normalised D. Melanogaster TFs', fontsize = 17)
plt.savefig('/Users/janadelmann/polybox/LabRotation1/figures/mean_info_content_normalised.png')
plt.show()
