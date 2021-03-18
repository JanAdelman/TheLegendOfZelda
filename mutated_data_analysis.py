import pandas as pd 
import numpy as np 
import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr
import pybedtools
from pybedtools import BedTool
import re

def mean(df):
    len = df.shape[1]

    count = 0
    for i in df:
        count += int(i)
    
    return np.round(count/len)

path = str(pathlib.Path(__file__).parent.absolute())

data = pd.read_csv(path + '/position_score_mutaed_sequences.csv', sep='\t')
original_zelda_count = data.iloc[:,1]
mutated_zelda_count = data.iloc[:,4]

mutation_occurances = np.where(data.iloc[:, 4] > data.iloc[:,1])
deleterious_occurances = np.where(data.iloc[:, 4] < data.iloc[:,1])

print(len(mutation_occurances[0]))
print(len(deleterious_occurances[0]))

# biasaway results:

# mononucleotide:

mononucleotide = pd.read_csv(path + '/BiasAwayResults/PWM_Scan/Mononucleotide/mononucleotide_count.csv', sep=',')
mean_mono = mean(mononucleotide)

ls = []
for i in mononucleotide:
    ls.append(int(i))
    
std_mono = np.std(ls)
print(std_mono)
print(mean_mono)

# dinucleotide

dinucleotide = pd.read_csv(path + '/BiasAwayResults/PWM_Scan/Dinucleotide/dinucleotide_count.csv', sep=',')
mean_di = mean(dinucleotide)

ls = []
for i in dinucleotide:
    ls.append(int(i))

std_di = np.std(ls)
print(std_di)
print(mean_di)

# random dna 
rand_dna = pd.read_csv(path + '/BiasAwayResults/random_dna_zelda_count.txt', sep=',')
mean_rand = mean(rand_dna)

ls = []
for i in rand_dna:
    ls.append(int(i))

std_rand = np.std(ls)
print(std_rand)
print(mean_rand)