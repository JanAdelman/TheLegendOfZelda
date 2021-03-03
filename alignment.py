import re 
import Bio
from Bio import SeqIO
import pandas as pd 
import glob

# In the following script the build ID47 of the d. Melanogaster genome was used 
# https://www.ncbi.nlm.nih.gov/genome?term=vih&cmd=DetailsSearch
# The following human genome assembly was used:
# https://www.ncbi.nlm.nih.gov/genome/?term=Homo+sapiens

path = "/Users/janadelmann/polybox/LabRotation1/drosophila_genome_split_files/*.fasta"

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
def get_count(file, zelda_motif):

    matches_dataframe = pd.DataFrame(columns=['sequence_id', 'start_position'])
    # read in the whole file at once 
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

    return matches_dataframe

zelda_motif = "CAGGTAG"

# Calculate matches for Drosophila:

drosophila_matches = get_count('/Users/janadelmann/polybox/LabRotation1/drosophila_genome.fa', zelda_motif)
print(drosophila_matches.shape)
drosophila_matches.to_csv('/Users/janadelmann/polybox/LabRotation1/matches_drosophila_genome', index=False)

# Calculate matches for H. Sapiens:

sapiens_matches = get_count('/Users/janadelmann/polybox/LabRotation1/Human_genome.fasta', zelda_motif)
print(sapiens_matches)
sapiens_matches.to_csv('/Users/janadelmann/polybox/LabRotation1/matches_human_genome', index=False)