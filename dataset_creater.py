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
        
