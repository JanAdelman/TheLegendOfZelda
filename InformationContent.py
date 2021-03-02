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

# https://www.ncbi.nlm.nih.gov/genome?term=vih&cmd=DetailsSearch
# median GC-content of ID-47, 42%

# The layout of the dataframe is: A C G T 
nucleo_distr = [0.29, 0.21, 0.21, 0.29]
nucleo_unif = 0.25

#==============================================================================#
# Calculate the Information content for each TF
#==============================================================================#

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

        # append the TF and the Infocontent to a data frame 
        
        df = pd.DataFrame([[name, info_tf]], columns = ['TF_Name', "Info"])
        TF_info = TF_info.append(df, ignore_index=True)
        TF_info = TF_info.sort_values(by='TF_Name')
        TF_info.to_csv('/Users/janadelmann/polybox/LabRotation1/information_content_TF', index=False)

#==============================================================================#
# Calculate basic statistics from the data TF
#==============================================================================#

# get the mean over all epxeriments 
mean_data = TF_info.groupby(by="TF_Name").mean()

# get the experiment with min informataion
min_data = TF_info.groupby(by="TF_Name").min()

# get the experiment with max information 
max_data = TF_info.groupby(by="TF_Name").max()


# Calculate average over all TF
mean = TF_info['Info'].mean()
print(mean)

# get information value for Zelda transcription factor 
mean_vfl = mean_data['Info']['vfl']
max_vfl = max_data['Info']['vfl']
min_vfl = min_data['Info']['vfl']

# get information value for GrainyHead transcription factor 
mean_grh = mean_data['Info']['grh']
max_grh = max_data['Info']['grh']
min_grh = min_data['Info']['grh']

#==============================================================================#
# Histogram Plot for different dataset compositions
#==============================================================================#

'''
# big plot with all criterions included 
colour = ['r', 'b', 'g']
#sns.histplot([mean_data['Info'], max_data['Info'], min_data['Info']], color=['r','b', 'g'], alpha=0.5)
data = [mean_data['Info'], max_data['Info'], min_data['Info'], ]
for d in data:
    sns.histplot(d)
plt.axvline(mean_vfl, color='darkblue')
plt.axvline(mean, color = 'red')
#sns.histplot(mean_data['Info'])
#sns.histplot(min_data['Info'])
#sns.histplot(max_data['Info'])
plt.xlabel('Information content (bits)', fontsize = 14)
plt.title('Information content of D. Melanogaster TFs')
plt.show()
'''
#fig, axs = plt.subplots(ncols=3)

# Plot of only average information 
sns.histplot(mean_data['Info'], alpha = 0.5)
plt.axvline(mean_vfl, color='darkblue')
plt.axvline(max_grh, color='green', alpha = 1)
plt.annotate('Zelda Information content', xy=(float(mean_vfl), 15), xytext=(float(mean_vfl)+2, 40), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Grainy Head Information content', xy=(float(mean_grh), 15), xytext=(float(mean_grh)+4.3, 55), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.ylabel("Count",fontsize=14)
plt.xlabel('Information content (bits)', fontsize = 14)
plt.title('Mean information content of D. Melanogaster TFs', fontsize = 17)
plt.savefig('/Users/janadelmann/polybox/LabRotation1/figures/mean_info_content.png')
plt.show()

# Plot of min information 
sns.histplot(min_data['Info'], alpha = 0.5)
plt.axvline(min_vfl, color='darkblue')
plt.axvline(min_grh, color='green')
plt.annotate('Zelda Information content', xy=(float(min_vfl), 15), xytext=(float(min_vfl)+2, 40), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Grainy Head Information content', xy=(float(min_grh), 15), xytext=(float(min_grh)+4.3, 55), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.xlabel('Information content (bits)', fontsize=14)
plt.title('Min. information content of D. Melanogaster TFs', fontsize=17)
plt.savefig('/Users/janadelmann/polybox/LabRotation1/figures/min_info_content.png')
plt.show()

# Plot of max information 
sns.histplot(max_data['Info'], alpha = 0.5)
plt.axvline(max_vfl, color='darkblue')
plt.axvline(max_grh, color='blue')
plt.annotate('Zelda Information content', xy=(float(max_vfl), 15), xytext=(float(max_vfl)+4, 40), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.annotate('Grainy Head Information content', xy=(float(min_grh), 15), xytext=(float(min_grh)+6, 55), arrowprops=dict(arrowstyle="->",facecolor='black'))
plt.xlabel('Information content (bits)', fontsize=14)
plt.title('Max. information content of D. Melanogaster TFs')
plt.savefig('/Users/janadelmann/polybox/LabRotation1/figures/max_info_content.png')
plt.show()
