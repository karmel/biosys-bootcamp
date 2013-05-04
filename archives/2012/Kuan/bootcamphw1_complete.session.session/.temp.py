import numpy as np
import scipy as sp
from pandas import DataFrame
from pandas.io import parsers
import os

os.chdir('/Users/kding/Desktop/bootcamp')

data = parsers.read_csv('data_set_HL60_U937_NB4_Jurkat.txt', sep='\t', header=1, index_col=1)