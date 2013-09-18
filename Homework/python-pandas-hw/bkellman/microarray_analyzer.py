import pandas,scipy,itertools,numpy,heapq,math,sets

#Input and parsing
data = pandas.read_csv('data_set_HL60_U937_NB4_Jurkat.txt','\t') 

a=pandas.Series(range(36))
a=a[a%2==0]
data2 = data[data.keys()[a]] 
hl60=data[data.keys()[a[1:5]]]
U937=data[data.keys()[a[5:9]]]
NB4=data[data.keys()[a[9:14]]]
Jurkat=data[data.keys()[a[14:18]]]
cellList = [hl60,U937,NB4,Jurkat]
cellnames = ['hl60','U937','NB4','Jurkat']

#A: gene count
print "GENE COUNT"
print len(sets.Set(data['Gene Description']))

#B: highest correlated combinations
c = lambda x,y: numpy.corrcoef(x,y)

highCors = list()
highCombs = list()
#iterate over cell types
for cell in cellList:
	highCor=0
	highComb=list()
	#iterate over all combinations of time points within cell types
	for subset in itertools.combinations(cell.keys(), 2):
		cor = abs(c(cell[subset[0]],cell[subset[1]]) [0,1])
		if cor > highCor:
			highCor = cor
			highComb = subset
	highCombs.append(highComb)
	highCors.append(highCor)
print "QUESTION B: HIGHEST CORRELATED TIME POINTS"
print highCombs
print highCors

#C: most similar cell lines
timeList = list()
#iterate by time points
for time in range(4):
	diffList = list()
	#within time points compair all pairs of cells
	for subset in itertools.combinations(cellList, 2):
		diff = abs(subset[0][subset[0].keys()[time]] - subset[1][subset[1].keys()[time]])
		diffList.append(diff)
	timeList.append(diffList)
		
#initialize combination names
combn=list()
for subset in itertools.combinations(cellnames, 2):
	combn.append('-'.join(subset))
	
#sum within cell comparisons accross time points
diffsum = list()
for i in range(6):
	print i
	diffsum.append( sum(timeList[0][i] + timeList[1][i] + timeList[2][i] + timeList[3][i]))

df=pandas.DataFrame({'names':pandas.Series(combn),'diffSums':pandas.Series(diffsum)})
print 'QUESTION C: MOST SIMILAR CELL TYPES'
print df.sort()

#iterate over cell types
#diffList = list()
#for subset in itertools.combinations(cellList, 2):
#	diff0 = 
#solution in text

#D:
#iterate through rows
sdList = list()
genes = list()
access = list()
count=0
for row in data2.iterrows():
	std = row[1][1::].std()
	sdList.append(std)
	genes.append(row[1][0])
	access.append(data.ix[count][1])
	count+=1

sdGenes=pandas.Series(sdList,index=genes)
sdAcc=pandas.Series(sdList,index=access)
#df = pandas.DataFrame({'genes':sdGenes,'access':pandas.Series(access)})
pandas.Series.sort(sdGenes)
pandas.Series.sort(sdAcc)
print "QUESTION D: MOST CONSISTANT EXPRESSIONS ACCROSS ALL CELL TYPES AND PERTERBATIONS"
print sdGenes[0:20]
print sdAcc[0:20]

#E: x>2-fold increase in expression 0->24hr
change = lambda new,old,old2: ((new-old)*numpy.sign(1.0,old2))/old
#iterate over cells
for cell in cellList:
	#shift data
	#cell['fold'] = cell[cell.keys()[3]] / cell[cell.keys()[0]]
	cell['fold'] = change( cell[cell.keys()[3]] , cell[cell.keys()[0]] , cell[cell.keys()[0]].copy() )
	cell['Access'] = data[data.keys()[1]]
	cell['Genes'] = data[data.keys()[0]]
	inc = cell[cell['fold']>2]
	inc = inc[inc['fold']!=float('Inf')]
	print 'QUESTION E: 2-FOLD INCREASE IN GENES FROM 0->24'
	print inc[[cell.keys()[0],cell.keys()[3],'Access','fold','Genes']].head()
	print inc['Access']

#F:differential expression between hl60 and u937
cell2 = pandas.DataFrame({ 'fold': change( hl60[hl60.keys()[0]] ,  U937[U937.keys()[0]] , U937[U937.keys()[0]].copy() ) ,
	'Access': data[data.keys()[1]],
	'Genes': data[data.keys()[0]] })

inc = cell2[cell2['fold']>2]
inc = inc[inc['fold']!=float('Inf')]
dec = cell2[cell2['fold']<-2]
dec = dec[dec['fold']!=float('Inf')]
print 'QUESTION F: DIFFERENTIAL EXPRESSION BETWEEN HL60 AND U937'
print 'increasing genes U937 -> HL60'
#print inc[['Access','fold','Genes']].head()
print inc['Access'].values[0:674]
print inc['Access'].values[674::]
#print inc['Access'].to_csv(sys.stdout)
print 'decreasing genes U937 -> HL60'
#print dec[['Access','fold','Genes']].head()
print dec['Access'].values
