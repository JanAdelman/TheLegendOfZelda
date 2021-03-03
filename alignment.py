import re 
import Bio
from Bio import SeqIO
import pandas as pd 
import glob

path = "/Users/janadelmann/polybox/LabRotation1/drosophila_genome_split_files/*.fasta"
zelda_motif = "CAGGTAG"

matches_dataframe = pd.DataFrame(columns=['sequence_id', 'start_position'])

#==========================================================#
# Find matches of the Zelda motif in the  Drosophila genome
#==========================================================#

# This function would read the split files separately 
'''
# read in all the fasta files 
for file in glob.glob(path):
    # get sequence and ID 
    for seq_record in SeqIO.parse(file, "fasta"):
        
        # Transform all entries to uppercase, transform to string for finditer function 
        seq_record.seq = str(seq_record.seq.upper())

        # find the pattern in the sequence 
        matches = re.finditer(zelda_motif, seq_record.seq)
        for m in matches: 
            base = m.group() 
            pos  = m.start() 

            # record all the found entries in a data frame 
            df = pd.DataFrame([[seq_record.id,m.start()]], columns=['sequence_id', 'start_position'])
            matches_dataframe = matches_dataframe.append(df,ignore_index=True)
'''

# read in the whole file at once 
for seq_record in SeqIO.parse('/Users/janadelmann/polybox/LabRotation1/drosophila_genome.fa', "fasta"):
        
        # Transform all entries to uppercase, transform to string for finditer function 
        seq_record.seq = str(seq_record.seq.upper())

        # find the pattern in the sequence 
        matches = re.finditer(zelda_motif, seq_record.seq)
        for m in matches: 
            base = m.group() 
            pos  = m.start() 

            # record all the found entries in a data frame 
            df = pd.DataFrame([[seq_record.id,m.start()]], columns=['sequence_id', 'start_position'])
            matches_dataframe = matches_dataframe.append(df,ignore_index=True)

print(matches_dataframe.shape)
matches_dataframe.to_csv('/Users/janadelmann/polybox/LabRotation1/matches_drosophila_genome', index=False)