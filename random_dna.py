import random
import pandas as pd
import sys
import glob
import os 
import pathlib

path = str(pathlib.Path(__file__).parent.absolute())

def random_sequennce(len):
    dna = ['A', 'C', 'G', 'T']
    sequence = ''
    for i in range(0, len):
        sequence += random.choice(dna)
    
    return sequence 

chromosome_dataset = pd.read_csv(path + '/chromFa/chrom_sizes_dm3.txt', sep='\t')
len = chromosome_dataset.iloc[:,1].sum()

sequence = random_sequennce(len)
sequence_name = 'random_dna_sequence'

ofile = open(path + "/random_dna_dm3.fa", "w")
ofile.write(">" + sequence_name + "\n" + sequence + "\n")
ofile.close()

