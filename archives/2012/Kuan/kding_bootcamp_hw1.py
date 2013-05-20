#Kuan-Fu Ding, Bootcamp - Python HW 1

'''I first import a set of packages that include functions used in script
(apparetly during cleaning/debugging I removed some of the functions and no longer used
some of the packages'''
# Yeah, be careful with that. If you have a tendency to accidentally break your own code,
# I highly, highly suggest you version all code using git or an equivalent.
import numpy as np
import scipy
from pandas import *
from pandas.io import *
import os
from scipy import *

# Where is my object?

'''change current working directory and load tab delimited file
		-set rownames to be Gene Accession IDs due to redundant genes
		-set the column header to first row
		-skip second and third blank rows
		-delete last row of NA values
'''

os.chdir('/Users/kding/Desktop/bootcamp')
data = parsers.read_csv('data_set_HL60_U937_NB4_Jurkat.txt', sep="\t", index_col=1, header=0, skiprows=[1,2], nrows = 7229)


'''To check unique genes, first extract Gene Discription columns and use the set function.
This lets us remove duplicated genes.  I then created a list of the duplicated genes and used
the len function to calculate the length'''
# Love the detailed comments, but I advise getting used to writing in passive sentences,
# or in sentences where the code itself is the actor.
# Bad for English, good for scientists.
unique = list(set(data['Gene Description']))
len(unique) ##5001

'''Extract data columns only - create new array with only expression data
(deleted calls and gene descriptions.)
Create subset of values for each cell type
 - I thought this was useful for within cell type coparisons
'''
values = data[data.columns[range(1,data.shape[1],2)]]
hl60 = values[values.columns[range(0,4)]]
u937 = values[values.columns[range(4,8)]]
nb4 = values[values.columns[range(8,13)]]
jurkat = values[values.columns[range(13,17)]]

'''Generate pairwise correlation matrix with corr function'''
hl60.corr()
# 0 and 4 hours have highest correlation
u937.corr()
# 0 and .5 hours have highest correlation
nb4.corr()
# 0 and 5.5 hours
jurkat.corr()
# 0 and 4 hours

'''Create array with expression data from 0 and 24 hours (these were the only common time pts)
Generate pairwise correlation matrix with corr function'''
comparecell = values[values.columns[[0,4,8,13]]]
comparecell.corr()  #u937 and jurkat are most similar at 0 hour
comparecellslater = values[values.columns[[3,7,10,16]]]
comparecellslater.corr()  #u937 and jurkat are most similar at 24 hour

'''Look at variance (or stdev) and sort from smallest to greatest.
The ten gene accession ids with the smallest variation across all samples are listed''' 
allvar = values.var(axis = 1)
allvar.sort
'''
10 smallest variance are:
H85478                   2.735294
H88261                   2.882353
M95586_r_i               3.632353
X73424_i                 3.779412
R94219                   4.183824
U39231                   4.235294
H82137                   4.264706
L20814                   4.279412
R39006                   4.360294
X86401                   4.404412
'''

'''I checked for two fold expression in two cases. The first was within cell type
only.  The second was within all four cell types.  The gene accession IDs an fold
change are listed.  Note: negative fold change is represented as decimal.  Take
the negative reciprocal for neg FC value.
-The first step in both cases was to check for present call.  Remove absent
probesets.  For within cell type, check for presence only in specified cell type.
For within all cell type, check for Present in all.
-Then, Generate fold change values and use logial operators to filter FC > abs(2)
'''
# Clever of you to figure out what those call columns were for. I didn't even consider them.
hl60present = data[(data[data.columns[2]] == 'P') & (data[data.columns[8]] == 'P')]
diffhl60 = hl60present[hl60present.columns[7]]/hl60present[hl60present.columns[1]]
diffhl60[np.logical_or((diffhl60>2), (diffhl60<.5))]   #704

u937present = data[(data[data.columns[10]] == 'P') & (data[data.columns[16]] == 'P')]
diffu937 = u937present[u937present.columns[15]]/u937present[u937present.columns[9]]
diffu937[np.logical_or((diffu937>2), (diffu937<.5))]   #511

nb4present = data[(data[data.columns[18]] == 'P') & (data[data.columns[22]] == 'P')]
diffnb4 = nb4present[nb4present.columns[21]]/nb4present[nb4present.columns[17]]
diffnb4[np.logical_or((diffnb4>2), (diffnb4<.5))]   #427


jurkatpresent = data[(data[data.columns[28]] == 'P') & (data[data.columns[34]] == 'P')]
diffjurkat = jurkatpresent[jurkatpresent.columns[33]]/jurkatpresent[jurkatpresent.columns[27]]
diffjurkat[np.logical_or((diffjurkat>2), (diffjurkat<.5))]   #281

presentall = data[(data[data.columns[2]] == 'P') & (data[data.columns[8]] == 'P') & (data[data.columns[10]] == 'P') & (data[data.columns[16]] == 'P')
& (data[data.columns[18]] == 'P') & (data[data.columns[22]] == 'P') & (data[data.columns[28]] == 'P') & (data[data.columns[34]] == 'P')]
diffhl60_all = presentall[presentall.columns[7]]/presentall[presentall.columns[1]]
diffu937_all = presentall[presentall.columns[15]]/presentall[presentall.columns[9]]
diffnb4_all = presentall[presentall.columns[21]]/presentall[presentall.columns[17]]
diffjurkat_all = presentall[presentall.columns[33]]/presentall[presentall.columns[27]]
diffhl60_all[np.logical_and(np.logical_and((np.logical_or((diffhl60_all>2), (diffhl60_all<.5))), (np.logical_or((diffu937_all>2), (diffu937_all<.5)))),
																													(np.logical_and((np.logical_or((diffnb4_all>2), (diffnb4_all<.5))),(np.logical_or((diffjurkat_all>2), (diffjurkat_all<.5))))))]

"""
Genes with greater than 2 fold change in all four cell types
Gene Accession Number
H42477                   3.142857
L06132                   0.309038
L20859                   3.805970
M15395                   4.689189
M37033                   8.028986
M61832                   0.289474
M88279                   0.138158
R56440                   0.358382
R61359                   3.136364
R61502                   0.343407
T70920_f                 0.236080
T70920_r_i               0.145973
U30498                   5.729032
X06956                   0.317308
"""

'''Checking for fold change at 0 hrs between hl60 and u937.  Similar procedure 
as above.
'''
hlupresent = data[(data[data.columns[2]] == 'P') & (data[data.columns[10]] == 'P')]
diffhlu = hlupresent[hlupresent.columns[9]]/hlupresent[hlupresent.columns[1]]
diffhluf = diffhlu[np.logical_or((diffhlu>2), (diffhlu<.5))]     ##913 accession IDS
diffhluf.to_csv('finalout.csv')