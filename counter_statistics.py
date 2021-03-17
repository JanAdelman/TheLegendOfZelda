import pandas as pd 
import glob
import os 
import pathlib
import pyranges as pr
import numpy as np 
import seaborn as sns 
import matplotlib
import matplotlib.pyplot as plt
from numpy import genfromtxt
import seaborn as sns

path = str(pathlib.Path(__file__).parent.absolute())
counts = pd.read_csv(path+'/counts.txt', delimiter=',').set_index('identity')
counts = counts.iloc[:,0:16]

counts.sum(axis = 0, skipna = True).plot.bar()
plt.title('Dinucleotide Counts D. Melanogaster')
plt.ylabel('counts')
plt.xlabel('Dinucleotides')
plt.show()

