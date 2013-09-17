#########################################################################################
#                                     IMPORTS                                           #
#########################################################################################

# Import division from the future

from __future__ import division

# Import argv

from sys import argv

# Import NumPy

import numpy as np

# Import pandas

import pandas
from pandas import Series
from pandas import DataFrame


#########################################################################################
#                                    DATA FRAME                                         #
#########################################################################################

# The solution takes the dataset as its argument, so we assign the argument to the 
# variable "dataset." We convert the input argument to a list so that lines are indexed.

dataset = list(open(argv[1], 'r'))

# The first and last lines of the input file contain unwanted HTML code, so we remove 
# them.

dataset.pop(0)
dataset.pop(len(dataset)-1)

# We also want to separate the datapoints on each line. They are separated by tab
# characters.

for line in range(len(dataset)):
    dataset[line] = dataset[line].replace(' call', '\tcall').split('\t')

# The resulting list of lists is converted easily to a NumPy array. Since the entries in
# the second column are unique while those in the first column are not, we would prefer 
# to exchange the two columns and use the Accession numbers to identify rows in the data 
# frame. We also do away with the unneeded "call" columns, convert strings to floating
# points where appropriate, and clean up one of the column titles.

dataset = np.array(dataset)
dataset[:,[0,1]] = dataset[:,[1,0]]
dataset = np.append( dataset[:,0:2], dataset[:,2:dataset.shape[1]-1:2], axis=1)
dataset[0,1] = 'Gene Description'

dataset[1:,2:] = dataset[1:,2:].astype(np.float)

# Now we can produce a dataframe.

dataframe = DataFrame(dataset[1:,1:],
                    index=dataset[1:,0],
                    columns=dataset[0,1:])

#########################################################################################
#                              COUNTING UNIQUE GENES                                    #
#########################################################################################

# Counting the number of unique genes is now easy.

genecount = len(set(dataframe['Gene Description']))

#########################################################################################
#                              TIME POINT CORRELATION                                   #
#########################################################################################

# We want to find out, for each cell type, which two time points are most tightly
# correlated.

# First, let's find out how many cell types there are and what they are called.

columns = dataframe.columns
celltypes = list(columns)

for column in range(len(columns)):
    celltypes[column] = dataframe.columns[column].split('_')[0]
    
celltypes.pop(0)
celltypes = list(set(celltypes))

