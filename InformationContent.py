# import necessary packages
import numpy as np
import os
import csv
import re 
import glob
from pathlib import Path
import pandas as pd

# get working directory
cwd = os.getcwd()

# create new directory for the TF files to be put in 
Path("/Users/janadelmann/polybox/LabRotation1/TF_PWM").mkdir(parents=True, exist_ok=True)

# read in the data from full csv file and add the each TF to a seaparate csv file
#with open('/Users/janadelmann/polybox/LabRotation1/B1H_TF_DATA/dummyfile.txt', newline='') as fullfile:
with open('/Users/janadelmann/polybox/LabRotation1/B1H_TF_DATA/PWM.txt', newline='') as fullfile:
    
    lines = fullfile.readlines()
    i = 1
    check = 0

    # create csv file for one TF 
    tf_csv = open(f"/Users/janadelmann/polybox/LabRotation1/TF_PWM/TF_{i}.csv", "w")

    # loop trough big data frame 
    for line in lines:

        # if line is empty one loop trough one of the TFs is finished
        # variable check assures that only at the second empty line a new file is created
        # in the data structure above after every TF two empyt lines exist 
        if line.strip() == '':
            if check == 1:
                i+=1
                tf_csv =  open(f"/Users/janadelmann/polybox/LabRotation1/TF_PWM/TF_{i}.csv", "w")
                check = 0
            
            else:
                check += 1
        
        # append TF data to csv file 
        else:
            tf_csv.write(line)  
        
            

def get_TF_data(TF_data):
    '''
    Function to return a pandas data frame for one TF, the data has not the same column number 
    for each row thus column are appended such that the dimensions match 
    '''

    with open(fname, newline='') as temp_f:
        # get No of columns in each line
        col_count = [ len(l.split("\t")) for l in temp_f.readlines() ]
           
        
    column_names = [i for i in range(0, max(col_count))]
    df = pd.read_csv(fname, header=None, delimiter="\t", names=column_names)
    name = str(df[1][1])
    print(name)
    df = df.drop([0,1,2,3,4,5,6])

    return df, name

path = "/Users/janadelmann/polybox/LabRotation1/TF_PWM/*.csv"
# get data for each TF
TF_info = pd.DataFrame(columns=['TF_Name', "Info"])

for fname in glob.glob(path):

    filesize = os.path.getsize(fname)
    
    # if file is empty, skip 
    if filesize == 0:
        pass
    else: 
        TF_data, name = get_TF_data(fname)

        print(TF_data)
        
        info_tf = 1 
        df = pd.DataFrame([[name, info_tf]], columns = ['TF_Name', "Info"])
        TF_info = TF_info.append(df, ignore_index=True)
        print(TF_info)

    '''
    with open(fname, newline='') as TF_file:
        
        lines = TF_file.readlines()
        
        p = re.compile('Gene')
        q = re.compile('Family')
        l = re.compile('TF Name')
        m = re.compile('Species')
        n = re.compile('Motif')
        o = re.compile('')

    

    for line in lines:

        if(p.search(line) or q.search(line) or l.search(line) or m.search(line) or n.search(line) or line.strip() == ''):
            continue   
        else:
            print(line)

    '''

