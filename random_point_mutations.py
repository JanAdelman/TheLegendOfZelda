import random
import numpy as np 
import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr
import pybedtools
from pybedtools import BedTool
import re
from random import choices
import Bio
from Bio import motifs

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

#import count_dinuclotides

path = str(pathlib.Path(__file__).parent.absolute())

# import the position weight matrix 
# pwm_zelda = np.loadtxt(path + '/vfl/vfl_motifs/vfl_full.txt', delimiter='\t')
# pwm_zelda = np.transpose(pwm_zelda)

with open(path + '/vfl/vfl_motifs/vfl_full.pfm') as handle:
    pwm_zelda = motifs.read(handle, "pfm")

# add pseudo count to avoid having 0 
pwm_zelda = pwm_zelda.counts.normalize(0.01)

# if needed print consensous sequence 
cons = pwm_zelda.degenerate_consensus

# background frequency of drosophila genome 
background = {'A':0.29,'C':0.21,'G':0.21,'T':0.29}

# calculate pssm (log odds)
pssm_zelda = pwm_zelda.log_odds(background=background)
mean = pssm_zelda.mean(background=background)
std = pssm_zelda.std(background=background)

# A uniform background is used if background is not specified. The mean is particularly important, 
# as its value is equal to the Kullback-Leibler divergence or relative entropy, and is a measure for the information 
# content of the motif compared to the background. As in Biopython the base-2 logarithm is used 
# in the calculation of the log-odds scores, the information content has units of bits.

#print("mean = %0.2f, standard deviation = %0.2f" % (mean, std))

# Selecting a threshold 
distribution = pssm_zelda.distribution(background=background, precision=10**4)

# calculate threshold (only calculate once, takes time)
#threshold = distribution.threshold_fpr(0.0001)
#threshold = 9.74
threshold = 9.74

# threshold = distribution.threshold_patser()
# print("%5.3f" % threshold)

def sample_row(full_df):
    '''
    This function returns a random row of a dataframe 
    Input:  pandas dataframe 
    '''

    return full_df.sample(n=1,axis='rows',replace=True)


def point_mutation(seq, num_mutations):
    '''
    This function performs random point mutations
    Input:  seq: The sequence to be mutated
            num_mutations: The number of point mutations that will occur       
    '''

    seq_len  = len(sequence)

    for repeats in range(num_mutations):

        # get a random position in the sequence
        rand_pos = random.randint(0, seq_len-1)
        # define a background distribution for the nucleotides 
        nucleotides = ['A', 'C', 'G', 'T']
        background_distr = [0.21, 0.29, 0.29, 0.21]
        
        # sample a nucleotide according to the background distribution 
        sampled_nucleotide = (str(random.choices(nucleotides, background_distr))).strip('[]').strip("''")
       
        # replace the nucleotide x at rand_pos by sampled_nucleotide
        seq = list(seq)
        seq[rand_pos] = sampled_nucleotide
        seq = ''.join(seq)

    return seq  

def pwm_scan(seq, pssm, threshold):
    '''
    This function performs a position weight matrix scan to find motifs in the sequence 
    '''
    position = []
    scores = []
    counter = 0
    for pos, score in pssm.search(seq, threshold=threshold):
        counter +=1
        position.append(pos) 
        scores.append(score)
        #print('positions mutated:  ', pos,'score mutated: ', score)
    
    return counter, position, scores 
  
# read in the dataframe 
colnames = ['Chrom', 'Start', 'End']
enhancer_info = pd.read_csv(path + '/tiles.csv')

#Index the drosophila genome by the row and return the sequence 
drosophila_genome = path + '/chromFa/dm3_all_chromosomes.fa'

sampled_df = pd.DataFrame()
# sample a random row from the df 
samples = 10000
# how many point mutations will be introduced 
num_mutations = 1

sampled_sequences = pd.DataFrame(columns=['original_seq', 'original_zelda_count', 'orig_pos', 'orig_score', 'mutated_seq', 'mutated_zelda_count', 
'mut_pos', 'mut_score'])

for repeats in range(samples): 

    # sample row to mutate 
    sampled_row = sample_row(enhancer_info)
    sampled_df = sampled_df.append(sampled_row)

    # convert panda df to bed file to be used to index fasta file and retrive the sequence 
    bed_format = sampled_row.iloc[:,1:4].to_numpy()

    np.savetxt(path + '/sample_for_mutation.bed', bed_format, delimiter='\t', fmt='%s')   
    get_sequence = pybedtools.BedTool(path + '/sample_for_mutation.bed')
    
    # get the sequence from the dm3 drosophila genome 
    sequence_to_mutate = get_sequence.sequence(fi=drosophila_genome)
    sequence = open(sequence_to_mutate.seqfn).read()
   
    # remove first line of fasta file 
    sequence = re.sub(r'>.*\n?', '', sequence, flags=re.MULTILINE)

    # test if in this sequence a zelda motif exists 

    #perform pwm scan with the sequence
    counter, position, scores = pwm_scan(sequence, pssm_zelda, threshold)

    # add point mutation to file 
    mutated_sequence = point_mutation(sequence, num_mutations)

    # scan the mutated sequence if zelda motif occured 
    mut_seq_motif_count, position_mut, scores_mut = pwm_scan(mutated_sequence, pssm_zelda, threshold)
 
    sampled_sequences = sampled_sequences.append({'original_seq': sequence, 'original_zelda_count': counter , 'mutated_seq': mutated_sequence, 
    'mutated_zelda_count': mut_seq_motif_count, 'orig_pos': position, 'orig_score': scores, 'mut_pos': position_mut, 'mut_score': scores_mut}, ignore_index=True)
    
    
sampled_sequences.to_csv(path + '/mutated_sequences.csv', sep='\t')

#print(sampled_sequences)      




