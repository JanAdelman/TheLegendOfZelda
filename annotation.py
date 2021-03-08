import re 
import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr

path = str(pathlib.Path(__file__).parent.absolute())
print(path)

# read in data
enhancer_info = pd.read_csv(str(path) + '/tiles.csv')
enhander_chromosome_names = enhancer_info.Chrosome.unique()

zelda_hits_full = pr.read_bed(str(path )+ '/TallyNumberOfSites/pwmscan_dm3_5002_7753.bed', as_df=True)
hits_chromosome_names = zelda_hits_full.Chromosome.unique()

chromosome_names = set.intersection(set(hits_chromosome_names), set(enhander_chromosome_names))

zelda_hits = zelda_hits_full['Start']

enhancer_zelda_motif_match = pd.DataFrame(columns = ['Chromosome', 'Enhander_start', "Enhander_end", "Zelda_site"])

for chromosome in chromosome_names:

    zelda_hits_per_chromosome = zelda_hits_full.loc[zelda_hits_full['Chromosome']==chromosome]
    enhancer_per_chromosome = enhancer_info.loc[enhancer_info['Chrosome']==chromosome]
    zelda_hits_per_chromosome = zelda_hits_per_chromosome['Start']

    for position in zelda_hits_per_chromosome:
        for i in range(0, enhancer_per_chromosome.shape[0]):
            if (position >= enhancer_per_chromosome.iloc[i][2] and position <= enhancer_per_chromosome.iloc[i][3]):

                df = pd.DataFrame([[chromosome, enhancer_per_chromosome.iloc[i][2], enhancer_per_chromosome.iloc[i][3], position]], columns = ['Chromosome', 'Enhander_start', "Enhander_end", "Zelda_site"])
                enhancer_zelda_motif_match = enhancer_zelda_motif_match.append(df, ignore_index=True)
                enhancer_zelda_motif_match.to_csv(path + '/enhander_zelda_matches', index=False)
                
                '''
                print('==========================')
                print('match found:')  
                print('at enhancer start position: ', enhander_per_chromosome.iloc[i][2])
                print('Zelda hit position: ', position)
                print('enhander end position: ', enhander_per_chromosome.iloc[i][2])
                print('on chromosome: ', chromosome)
                print('==========================')
                '''
       

