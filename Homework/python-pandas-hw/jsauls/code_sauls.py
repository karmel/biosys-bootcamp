# Get dataframe
from pandas import *

# Import data from text file
data = open('/Users/jt_temp/Desktop/data_set_HL60_U937_NB4_Jurkat.txt','r')
df = read_csv(data,index_col=0,sep='\t')

rows = df.index
rows = len(rows)

print 'A. There are',rows,'genes.'

# Important columns for easy grabbing
cols = ['Gene Accession Number','HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs',
        'U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs','NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs',
        'NB4_48_hrs','NB4_72_hrs','Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']

print 'Most similar timepoint by cell.'
# Import numpy
import numpy as np

# Function that finds Euclidian distance between gene arrays.
def gene_dist(arrayname1,arrayname2):
    arrayval1 = df[arrayname1]
    #print(str(type(arrayval1)))
    arrayval2 = df[arrayname2]
    eu_dist = np.linalg.norm(arrayval1-arrayval2)
    return eu_dist

# Ranges of the cell types in my above cols
cell_ranges = [[1,4],[5,8],[9,13],[14,17]]
for cell in cell_ranges:
    # Init dist
    dist_low = -1
    n = cell[0]
    m = cell[1]
    for start_array in range(n,m): # The range here is the timepoint for one cell. 
        for compare_array in range(n+1,m+1):
            # Call the function for every combination.
            dist_temp = gene_dist(cols[start_array],cols[compare_array])
            #print dist_temp,str(start_array),str(compare_array)
            # Record the first distance
            if dist_low == -1:
                start_index = start_array
                compare_index = compare_array
                dist_low = dist_temp
            # Record lower distances
            if dist_temp < dist_low:
                start_index = start_array
                compare_index = compare_array
                dist_low = dist_temp
        # Update compare array 
        n += 1
    print 'The most similar timepoints are between',cols[start_index],'and',cols[compare_index],'at',str(dist_low)

# Problem C
# I know what I want to do here: first compare the transcription profiles of each gene pairwise between the
# cells. This will result in 3 arrays (HL60 vs U937, HL60 vs Jurkat, and U937 vs Jurkat) of the norm distance 
# between the gene's transcription profiles. I would omit the other dataset because its timepoints are out of 
# sync and I'm not sure how to weight for that. Then I would just compute the sum of the abs value of those 
# and whichever one is the lowest, I would call those most similar. But I can't for the life of me get the 
# columns arranged right just to do the first calculation!

# Declare cells, skip NB4 because the timepoints ars so wacky
HL60 = df[[cols[1],cols[2],cols[3],cols[4]]]
HL60 = HL60.T
HL60 = HL60.iloc[0:4,0:7229]
print HL60

U937 = df[[cols[5],cols[6],cols[7],cols[8]]]
U937 = U937.T
U937 = U937.iloc[0:4,0:7229]
print U937

Jurkat = df[[cols[14],cols[15],cols[16],cols[17]]]
Jurkat = Jurkat.T
Jurkat = Jurkat.iloc[0:4,0:7229]
print Jurkat

# Comprare trancription profiles
def comp_prof(prof1,prof2):
    eu_dist = np.linalg.norm(arrayval1-arrayval2)
    return eu_dist

prof1 = HL60['BioB (spiked control)']
print prof1
prof2 = U937['BioC (spiked control)']
print prof2
prof1-prof2


print 'Least variable genes for calibration.'

# Data containing df
df_data = df[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs',
        'U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs','NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs',
        'NB4_48_hrs','NB4_72_hrs','Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']]

# Finds standard deviation across rows
std_array = df_data.std(axis=1)

# Minimum indexes, use loop that then drops the index to find the next.
for minimum in range(0,10):
    min_index = np.argmin(std_array)
    print str(min_index),df['Gene Accession Number'].ix[min_index]
    std_array.ix[min_index] = 1000 # Take in out of the running


print 'Genes with 2-fold higher expression at 24h.'

df_double = df[(df['HL60_24_hrs'].abs() >= df['HL60_0_hrs'].abs()*2)
               & (df['U937_24_hrs'].abs() >= df['U937_0_hrs'].abs()*2)
               & (df['NB4_24_hrs'].abs() >= df['NB4_0_hrs'].abs()*2)
               & (df['Jurkat_24_hrs'].abs() >= df['Jurkat_0_hrs'].abs()*2)]

df_double = df_double['Gene Accession Number']
print df_double


print 'Genes that are diffenentially regulated between HL60 and U937 at 0h'

df_double = df[((df['HL60_0_hrs'].abs() / df['U937_0_hrs'].abs()) > 2)
               | ((df['HL60_0_hrs'].abs() / df['U937_0_hrs'].abs()) > 2)]

df_double = df_double['Gene Accession Number']
print df_double

f = open('/Users/jt_temp/Desktop/output.txt', 'w')
f.write(df_double)
f.close()


