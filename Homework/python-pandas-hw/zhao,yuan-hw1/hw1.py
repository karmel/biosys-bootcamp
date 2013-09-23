import pandas
import itertools
import numpy
from numpy import std
from itertools import combinations
from pandas import Series
from pandas import DataFrame
import csv

out = open('accession.txt', "w")

df = pandas.read_csv('data_set_HL60_U937_NB4_Jurkat.csv')

#a
numgenes = len(set(df['Gene Description']))
print numgenes
print '\n'

#b
cells_b = (('HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs'),('U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs'),('NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs'),('Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs'))

all_diff = {}

for celltype in cells_b:
	combo_b = combinations(celltype,2)
	cell_diff = {}
	for comparison in combo_b:
		diff = df[comparison[0]]-df[comparison[1]]
		mean = abs(diff.mean())
		cell_diff[str(comparison[0])+' vs. '+str(comparison[1])] = mean
	min_b = min(cell_diff, key=cell_diff.get)
	all_diff[min_b]=cell_diff[min_b]
print all_diff
print '\n'

#c
cells_c = ('HL60_24_hrs','U937_24_hrs','NB4_24_hrs','Jurkat_24_hrs')
combo_c = combinations(cells_c,2)
most_similar = {}
for comparison in combo_c:
	diff = df[comparison[0]]-df[comparison[1]]
	mean = abs(diff.mean())
	most_similar[str(comparison[0])+' vs. '+str(comparison[1])] = mean
min_c = min(most_similar, key=most_similar.get)


print min_c
print most_similar[min_c]
print '\n'

#D

#numentries = len(set(df['Gene Accession Number']))
#
#standards = {}
#
#cell_indices = ((2,4,6,8),(10,12,14,16),(18,20,22,24,26),(28,30,32,34))
#for brosef in cell_indices:
#	for gene in range(0,numentries):
#		row = df.ix[gene]
#		stdval = ()
#		for time in brosef:
#			stdval += (row[time],)
#		stdev = numpy.std(stdval)
#		stdev = abs(stdev)
#		standards['gene: '+str(row[0])+'; cell column indices: '+str(brosef)]=stdev
#
#min_d = min(standards, key=standards.get)
#print min_d
#print standards[min_d]


numentries = len(set(df['Gene Accession Number']))

standards = []

cell_indices = ((2,4,6,8),(10,12,14,16),(18,20,22,24,26),(28,30,32,34))
for brosef in cell_indices:
	for gene in range(0,numentries):
		row = df.ix[gene]
		stdval = ()
		for time in brosef:
			stdval += (row[time],)
		stdev = numpy.std(stdval)
		stdev = abs(stdev)
		standards.append(('gene: '+str(row[0])+'; cell column indices: '+str(brosef),stdev))

standards = sorted(standards, key=lambda x: x[1])
standards = standards[:10]

print standards
print '\n'



#E



multiples = []

cell_indices = ((2,8),(10,16),(18,22),(28,34))
for brosef in cell_indices:
	for gene in range(0,numentries):
		row = df.ix[gene]
		multiple = abs(float(brosef[1])/float(brosef[0]))
		if multiple >= 2:
			multiples.append(str(row[0]))
		else:
			continue

multiples2 = []
for everybody in multiples:
	count = multiples.count(everybody)
	if count == 4:
		multiples2.append(everybody)
	else:
		continue

multiples2 = set(multiples2)
print multiples2
print len(multiples2)
print '\n'


#F

yo = []
for gene in range(0,numentries):
	row = df.ix[gene]
	if 2.*row[2] <= row[10]:
		yo.append(str(row[1]))
	elif 0.5*row[2] >= row[10]:
		yo.append(str(row[1]))
	else:
		continue

print len(yo)
for aasdfasdf in yo:
	out.writelines(str(aasdfasdf)+'\n')