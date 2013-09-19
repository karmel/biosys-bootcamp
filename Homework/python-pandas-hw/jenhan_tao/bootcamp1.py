import scipy
import pandas
import numpy
import sys

# read in data, the file is given as the first argument
file = sys.argv[1]
# read in file
data = pandas.read_csv(file,sep='\t')

# print name
print('Jenhan Tao\nBioinformatics Bootcamp Assignment 1\n')
print('Question A')

# determine the answer to Question A
answerA = len(set(data["Gene Description"]))
print('Based on the gene descriptions, there are ' + str(answerA) + ' genes')

# determine the answer to Question B
print('\nQuestion B')
HL6columns0 = ('HL60_0_hrs', 'HL60_0.5_hrs', 'HL60_4_hrs', 'HL60_24_hrs')
U937Columns = ('U937_0_hrs', 'U937_0.5_hrs', 'U937_4_hrs', 'U937_24_hrs')
NB4Columns = ('NB4_0_hrs', 'NB4_5.5_hrs', 'NB4_24_hrs', 'NB4_48_hrs', 'NB4_72_hrs')
JurkatColumns = ('Jurkat_0_hrs', 'Jurkat_0.5_hrs', 'Jurkat_4_hrs', 'Jurkat_24_hrs')
columns = (HL6columns0, U937Columns, NB4Columns, JurkatColumns)

for i in range(0,len(columns)):
    differencePairDict = {}
    for j in range(0,len(columns[i])):
        for k in range(j+1, len(columns[i])):
            difference = abs((data[columns[i][j]] - data[columns[i][k]]).mean())
            differencePairDict[difference] = columns[i][j] + " " +columns[i][k]
    minKey = min(differencePairDict.keys())
    print("For " +columns[i][0][:columns[i][0].index("_")] + " the most correlated time points are " + differencePairDict[minKey].split("_")[1] + " hrs and " + differencePairDict[minKey].split("_")[3]+" hrs")

# determine the answer to Question C
minPairValue = 1000;
minPair = "";
columns = ("HL60_24_hrs", "U937_24_hrs", "NB4_24_hrs", "Jurkat_24_hrs")
print('\nQuestion C')
for i in range(0,len(columns)):
    differencePairDict = {}
    for j in range(i+1,len(columns)):
         value = abs((data[columns[i]] - data[columns[j]]).mean())
         if value < minPairValue:
             minPairValue = value
             minPair = columns[i]+"_" + columns[j]
print("The most similar cell types are " + minPair.split("_")[0] + " and " + minPair.split("_")[3])

# determine the answer to Question D
print('\nQuestion D')
columns = (HL6columns0, U937Columns, NB4Columns, JurkatColumns)
sdDict = {}
avgDict = {}
for i in range(0, len(columns)):
    for j in range(0, len(data)):	
        values = [];
        for k in range(0, len(columns[i])):
	    values.append(data[columns[i][k]][j])
        sd = numpy.std(values)
        if not data["Gene Description"][j] in sdDict.keys():
            sdDict[data["Gene Description"][j]] = [];    
        sdDict[data["Gene Description"][j]].append(sd)
for key in sdDict:
    avgDict[numpy.mean(sdDict[key])] = key 

scores = avgDict.keys()
scores = sorted(scores)[:10]
genes =[]
for score in scores:
    genes.append(avgDict[score])
print("Ten good candiate genes are:\n"+str(genes))

# determine the answer to Question E
print('\nQuestion E')
columns24 = ('HL60_24_hrs', 'U937_24_hrs', 'NB4_24_hrs', 'Jurkat_24_hrs')
columns0 = ('HL60_0_hrs', 'U937_0_hrs', 'NB4_0_hrs', 'Jurkat_0_hrs')
winners = []
for i in range(0,len(data)):
    counter =0
    for j in range(0, len(columns24)):
        if data[columns24[j]][i] >= 2*data[columns0[j]][i]:
            counter = counter + 1
    if counter == 4:
        winners.append(data["Gene Description"][i])
print("The following genes satisfy the condition specified in the problem: \n"+str(set(winners)))

# determine the answer to Question F
print('\nQuestion F')
winners = []
for i in range(0,len(data)):
    counter =0
    if data["HL60_0_hrs"][i] >= 2*data["U937_0_hrs"][i] or data["HL60_0_hrs"][i] <= 1/2*data["U937_0_hrs"][i]:
        winners.append(data["Gene Accession Number"][i])
print("The following genes satisfy the condition specified in the problem: \n"+str(set(winners)))



# determine the answer to Question G
print('\nQuestion G')
# print out for copy and paste
#for winner in winners:
#    print(winner)
print('Yes there are enriched terms')
