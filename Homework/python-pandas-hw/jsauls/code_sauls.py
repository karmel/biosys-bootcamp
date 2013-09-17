# Import modules
from pandas import *
import numpy as np

# Import data from text file
data = open('/Users/jt_temp/Desktop/data_set_HL60_U937_NB4_Jurkat.txt','r')
df = read_csv(data,index_col=0,sep='\t')
rows = df.index
rows = len(rows)

### Problem A
# Reimporting because I named the index column for the df I use elsewhere.
data = open('/Users/jt_temp/Desktop/data_set_HL60_U937_NB4_Jurkat.txt','r')
dNames = read_csv(data,sep='\t')
uniq = dNames['Gene Description'].nunique() 
print 'A. There are',uniq,'genes.\n'

### Problem B
print 'B. Most similar timepoint by cell.'

# Important columns for easy grabbing
cols = ['Gene Accession Number','HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs',
        'U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs','NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs',
        'NB4_48_hrs','NB4_72_hrs','Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']

# Function that finds Euclidian distance between gene arrays.
def gene_dist(arrayname1,arrayname2):
    '''Takes two column names from the df and computes norm of the arrays.''' 
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
print '\nC.'
# I know what I want to do here: first compare the transcription profiles of each gene pairwise between the
# cells. This will result in 3 arrays (HL60 vs U937, HL60 vs Jurkat, and U937 vs Jurkat) of the norm distance 
# between the gene's transcription profiles. I would omit the other dataset because its timepoints are out of 
# sync and I'm not sure how to weight for that. Then I would just compute the sum of the abs value of those 
# and whichever one is the lowest, I would call those most similar.

# Define Cells
HL60 = df[[cols[1],cols[2],cols[3],cols[4]]]
U937 = df[[cols[5],cols[6],cols[7],cols[8]]]
Jurkat = df[[cols[14],cols[15],cols[16],cols[17]]]

# Function to comute norm between to transcripion profiles
def prof_dist(p1,p2):
    '''Takes two arrays (or annoying columns) and computes norm.'''
    sub_array = [0,0,0,0]
    # Subtract arrays... because I can't transpose them
    for i in range(0,4):
        sub_array[i] = p1[i] - p2[i]
    # Calculate distance.
    sub_array_dist = np.linalg.norm(sub_array)
    return sub_array_dist

# Calculate vector of diffencenes of transcription proflies by gene for 2 cell types.
gene_distances = np.zeros(rows)
for profile in range(0,rows):
    #print profile
    prof1 = HL60.ix[profile]
    prof2 = U937.ix[profile]
    gene_distances[profile] = prof_dist(prof1,prof2)

# Use sum of abs val of array as 'similarity score'
HL60vsU937 = np.max(np.absolute(gene_distances))
print 'Distance between HL60 and U936 is',HL60vsU937

# Repeat
# Calculate vector of diffencenes of transcription proflies by gene for 2 cell types.
gene_distances = np.zeros(rows)
for profile in range(0,rows):
    #print profile
    prof1 = HL60.ix[profile]
    prof2 = Jurkat.ix[profile]
    gene_distances[profile] = prof_dist(prof1,prof2)

# Use sum of abs val of array as 'similarity score'
HL60vsJurkat = np.max(np.absolute(gene_distances))
print 'Distance between HL60 and Jurkat is',HL60vsJurkat

# Calculate vector of diffencenes of transcription proflies by gene for 2 cell types.
gene_distances = np.zeros(rows)
for profile in range(0,rows):
    #print profile
    prof1 = U937.ix[profile]
    prof2 = Jurkat.ix[profile]
    gene_distances[profile] = prof_dist(prof1,prof2)

# Use sum of abs val of array as 'similarity score'
U937vsJurkat = np.max(np.absolute(gene_distances))
print 'Distance between U937 and Jurkat is',U937vsJurkat

### Problem D
print '\nD. Least variable genes for calibration.'

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

### Problem E
print '\nE. Genes with 2-fold higher expression at 24h.'

df_double = df[(df['HL60_24_hrs'].abs() >= df['HL60_0_hrs'].abs()*2)
               & (df['U937_24_hrs'].abs() >= df['U937_0_hrs'].abs()*2)
               & (df['NB4_24_hrs'].abs() >= df['NB4_0_hrs'].abs()*2)
               & (df['Jurkat_24_hrs'].abs() >= df['Jurkat_0_hrs'].abs()*2)]

df_double = df_double['Gene Accession Number']
print df_double

### Problem F
print '\nF. Genes that are diffenentially regulated between HL60 and U937 at 0h'

df_double = df[((df['HL60_0_hrs'].abs() / df['U937_0_hrs'].abs()) > 2)
               | ((df['HL60_0_hrs'].abs() / df['U937_0_hrs'].abs()) > 2)]

df_double = df_double['Gene Accession Number']
print df_double

f = open('/Users/jt_temp/Desktop/output.txt', 'w')
f.write(df_double)
f.close()


