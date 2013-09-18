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
#                                   SORTING DATA                                        #
#########################################################################################

# We want to find out, for each cell type, which two time points are most tightly
# correlated.

# First, let's find out how many cell types there are and what they are called. We create
# a list called "celltypes" which at first contains all of the column titles from our 
# dataframe. We shrink it to a list of the cell types used in this dataset.

columns = list(dataframe.columns)
celltypes = list(dataframe.columns)

for column in range(len(columns)):

    celltypes[column] = columns[column].split('_')[0]
    
celltypes.pop(0)
celltypes = list(set(celltypes))
celltypecount = len(celltypes)

# Now that we know all of the cell types, we can identify which columns refer to each cell
# type and sort the column titles into separate lists. We create a list of lists called 
# "celltype_sorting_bin" which will hold the separated columns.

# First, we initialize the list, ensuring it has the correct number of empt entries.

celltype_sorting_bin = []

for celltype_counter in range(celltypecount):

    celltype_sorting_bin.append([])
    
# Then we scan the dataframe column titles, identify which columns belong to which cell 
# type, and place them in the appropriate entry of the "sorting bin."

for celltypeindex in range(celltypecount):

    for column in columns:
    
        if column[:len(celltypes[celltypeindex])] == celltypes[celltypeindex]:
        
            celltype_sorting_bin[celltypeindex].append(column)


#########################################################################################
#                              COUNTING UNIQUE GENES                                    #
#########################################################################################

# Counting the number of unique genes is easy.

genecount = len(set(dataframe['Gene Description']))

#########################################################################################
#                              TIME POINT CORRELATION                                   #
#########################################################################################

# We want to find, for each cell type, the two time points for which the data is most
# closely correlated.

# First, we initialize a dictionary which will store the results.

TPC_result_dict = {}

# Since the number of cell types in the data is arbitrary, we start a for loop to operate
# over them.

for celltypeindex in range(celltypecount):

    # For each cell type, we need to find out how many time points have been recorded. 
    # We also initialize a matrix to store correlation coefficients.

    timepointcount = len(celltype_sorting_bin[celltypeindex])
    corr_value_storage = np.zeros((timepointcount, timepointcount))
    
    # We now want to find a correlation coefficient for each pair of distinct time points
    # and record it in the storage matrix, which requires a double for loop. Since we 
    # only want to check each pair once, we also add an if statement.
    
    for timepointindex in range(timepointcount):
    
        for timepointjndex in range(timepointcount):
        
            if timepointindex > timepointjndex:
            
                # We are now looking at a specific pair of time points. We can find their
                # correlation coefficient ("corr_value",) and the value of the
                # greatest coefficient found before the current pair ("max_corr.")
            
                corr_value = dataframe[celltype_sorting_bin[celltypeindex][timepointindex]].corr(dataframe[celltype_sorting_bin[celltypeindex][timepointjndex]])
                max_corr = np.amax(np.absolute(corr_value_storage))
                
                # Next, we check if the coefficient of the current time point pair greater
                # than the previous maximum. If so, we record it as the most highly 
                # correlated pair. We record it as a list containing a list, which makes
                # things more convenient if there is a tie.
                
                if corr_value > max_corr:
                    
                    most_correlated_times = [[celltype_sorting_bin[celltypeindex][timepointjndex], celltype_sorting_bin[celltypeindex][timepointindex]]]
                
                # In the case of a tie, we append the new pair so that both (or all) 
                # maximally correlated pairs are recorded.
                
                if corr_value == max_corr:
                    
                    most_correlated_times.append([celltype_sorting_bin[celltypeindex][timepointjndex], celltype_sorting_bin[celltypeindex][timepointindex]])
                
                # Lastly, we record the current correlation coefficient in the storage
                # matrix.
                
                corr_value_storage[timepointindex, timepointjndex] = corr_value
                
    TPC_result_dict[celltypes[celltypeindex]] = most_correlated_times
    
#########################################################################################
#                                 PRINTING RESULTS                                      #
#########################################################################################    

print('\nFor each cell type, the most highly correlated time points are:\n')

for key in TPC_result_dict:

    print('%s: %r' % (key, TPC_result_dict[key]))
    
print('\n')
