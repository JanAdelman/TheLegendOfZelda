import pandas as pd 
import glob
import os 
import pathlib
import numpy as np
import pybedtools
import pyfastx
from pybedtools import BedTool

path = str(pathlib.Path(__file__).parent.absolute())

enhancer_info = pd.read_csv(str(path) + '/tiles.csv')
matches = pd.read_csv(str(path) + '/zelda_enhancer_match_full_motif_0.0001')
enhancer_info = enhancer_info.iloc[:,1:4].to_numpy()

np.savetxt(path + '/enhancer_info.bed', enhancer_info, delimiter='\t', fmt='%s')   

# load the enhancer info data with the coordinates of the enhancer regions
# load the number of sites from the PWM scan 
enhancer_coordinates = pybedtools.BedTool('enhancer_info.bed')
zelda_site_coordinates = pybedtools.BedTool(path + '/TallyNumberOfSites/pwmscan_dm3_cutoff_prob_0.0001.bed')

# get the intersection of the the two files this will give all enhancer locations where zelda motif can be found 
intersection_zelda_sites = enhancer_coordinates.intersect(zelda_site_coordinates, wa=True)
# save the region in a bed file with the coordinates from the enhancer file 
intersection_zelda_sites.saveas(path + '/zelda_enhancer_match_bedtools_0.0001.bed')

# get the sequences from the enhancer regions. Use the coordinates to retrive the dna sequence  
intersection_zelda_sites = pybedtools.BedTool(path + '/zelda_enhancer_match_bedtools_0.0001.bed')
drosophila_genome = BedTool(path + '/chromFa/dm3_all_chromosomes.fa')
enhancer_sequences_with_zelda_motif = intersection_zelda_sites.sequence(fi=drosophila_genome)

sequence_of_enhancer_regions = enhancer_coordinates.sequence(fi=drosophila_genome)

# save the enhancer regions to use for dinucleotide count 
with open(path + '/enhancer_regions_dm3.fa', 'w') as enhancer_regions_dm3:
    enhancer_regions_dm3.write(open(sequence_of_enhancer_regions.seqfn).read())


# load the coordinates of the non enhancer regions then again retrive the sequences 
complement_enhancer_coordiantes = pybedtools.BedTool(path +'/complement_enhancer_regions.bed')
complement_sequence_of_enhancer_regions = complement_enhancer_coordiantes.sequence(fi=drosophila_genome)

with open(path + '/complement_enhancer_regions_dm3.fa', 'w') as complement_enhancer_regions_dm3:
    complement_enhancer_regions_dm3.write(open(sequence_of_enhancer_regions.seqfn).read())
