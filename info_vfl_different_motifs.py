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

path = "/Users/janadelmann/polybox/LabRotation1/vfl/vfl_motifs/*.txt"
nucleo_distr = [0.29, 0.21, 0.21, 0.29]
nucleo_unif = 0.25

# get data for each TF

TF_info = pd.DataFrame(columns=['TF_Name', "Info"])


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

num = 0 

for fname in glob.glob(path):
    
    filesize = os.path.getsize(fname)
    
    # if file is empty, skip 
    if filesize == 0:
        pass
    else: 
        # get the PVM Matrix + the name of the TF 
        pvm_matrix = np.loadtxt(fname, dtype=float, delimiter='\t')

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
    
    print(fname)
    print(info_tf)