import csv
import argparse
import collections
from Bio import SeqIO
import pandas as pd
import sys

def parseArgs():
        """Parse commandline arguments"""

        import argparse
        try:
                parser = argparse.ArgumentParser(description='Count dinucleotide frequency in a multifasta and write the output to a CSV.')
                parser.add_argument('--infile',
                                    action='store',
                                    help='The HHpred output file to parse.')
                parser.add_argument('--countsfile',
                                    action='store',
                                    help='Output file to store counts in.')
        except:
                print("An exception occured with argument parsing. Check your provided options.")
                traceback.print_exc()

        return parser.parse_args()

args = parseArgs() 

filename = sys.argv[4]

dinucleotides = ['AA','AT','AC','AG','TT','TA','TC','TG','CC','CA','CT','CG','GG','GA','GT','GC']
ratios = [('AA','TT'),('AC','GT'),('AG','CT'),('TT','AA'),('TC','GA'),('TG','CA'),('CC','GG'),('CA','TG'),('CT','AG'),('GG','CC'),('GA','TC'),('GT','AC')]

df = pd.DataFrame(columns = dinucleotides)

with open(args.infile, 'r') as ifh, open(args.countsfile, 'w') as cfh:
    
    for record in SeqIO.parse(ifh, 'fasta'):
        
        count = collections.defaultdict(int)
        id = record.id
        
        for i in range(len(record.seq)-2): 
            count[record.seq[i:i+2].upper()] += 1

        df_new  = pd.DataFrame([count.values()], columns=count.keys())
        df = pd.concat([df, df_new], axis =0)
        df = df.assign(identity = id) 
        df = df.set_index('identity')

df.to_csv(filename)

        