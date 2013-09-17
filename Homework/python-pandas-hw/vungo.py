import pandas
import numpy as np
from sys import argv

data = pandas.read_csv('data_set_HL60_U937_NB4_Jurkat.csv')
print 'added data table'

# PART (A) number of genes
data_nodup = data.drop_duplicates(cols='Gene Description', take_last=True)
data_nodup.shape[0]
print data_nodup.shape[0]

"""

#PART (B) which time points are most highly correlated
#CELL TYPE HL60

HL160_comparison_1v2 = np.corrcoef([data['HL60_0_hrs'],data['HL60_0.5_hrs']])
HL160_comparison_1v2corr = HL160_comparison_1v2[0,1]

HL160_comparison_1v3 = np.corrcoef([data['HL60_0_hrs'],data['HL60_4_hrs']])
HL160_comparison_1v3corr = HL160_comparison_1v3[0,1]

HL160_comparison_1v4 = np.corrcoef([data['HL60_0_hrs'],data['HL60_24_hrs']])
HL160_comparison_1v4corr = HL160_comparison_1v4[0,1]

HL160_comparison_2v3 = np.corrcoef([data['HL60_0.5_hrs'],data['HL60_4_hrs']])
HL160_comparison_2v3corr = HL160_comparison_2v3[0,1]

HL160_comparison_2v4 = np.corrcoef([data['HL60_0.5_hrs'],data['HL60_24_hrs']])
HL160_comparison_2v4corr = HL160_comparison_2v4[0,1]

HL160_comparison_3v4  = np.corrcoef([data['HL60_4_hrs'],data['HL60_24_hrs']])
HL160_comparison_3v4corr = HL160_comparison_3v4[0,1]

a = (('0v0.5', HL160_comparison_1v2corr),('0v4', HL160_comparison_1v3corr),('0v24', HL160_comparison_1v4corr),('0.5v4', HL160_comparison_2v3corr),('0.5v24', HL160_comparison_2v4corr),('4vs24', HL160_comparison_3v4corr))


print max(a)


# CELL TYPE Jurkat
jurkat_comparison_1v2 = np.corrcoef(data['Jurkat_0_hrs'],data['Jurkat_0.5_hrs'])
jurkat_comparison_1v3 = np.corrcoef(data['Jurkat_0_hrs'],data['Jurkat_4_hrs'])
jurkat_comparison_1v4 = np.corrcoef(data['Jurkat_0_hrs'],data['Jurkat_24_hrs call'])
jurkat_comparison_2v3 = np.corrcoef(data['Jurkat_0.5_hrs'],data['Jurkat_4_hrs'])
jurkat_comparison_2v4 = np.corrcoef(data['Jurkat_0.5_hrs'],data['Jurkat_24_hrs call'])
jurkat_comparison_3v4 = np.corrcoef(data['Jurkat_4_hrs'],data['Jurkat_24_hrs call'])

jurkat_coefficients = (('1v2',jurkat_comparison_1v2[0,1]),('1v3',jurkat_comparison_1v3[0,1]),('1v4',jurkat_comparison_1v4[0,1]),('2v3',jurkat_comparison_2v3[0,1]),('2v4',jurkat_comparison_2v4[0,1]),('3v4',jurkat_comparison_3v4[0,1]))

print 'for Jurkat' ,max(jurkat_coefficients)


#CELL TYPE Jurkat

U937_compar_1v2= np.corrcoef([data['U937_0_hrs'],data['U937_0.5_hrs']])
U937_comparison_1v2corr = U937_compar_1v2[0,1]

U937_compar_1v3=
np.corrcoef([data['U937_0_hrs'],data['U937_4_hrs']])
U937_comparison_1v3corr = U937_compar_1v3[0,1]


U937_compar_1v4=
np.corrcoef([data['U937_0_hrs'],data['U937_24_hrs']])
U937_comparison_1v4corr = U937_compar_1v4[0,1]

U937_compar_2v3=
np.corrcoef([data['U937_0.5_hrs'],data['U937_4_hrs']])
U937_comparison_2v3corr = U937_compar_2v3[0,1]

U937_compar_2v4=
np.corrcoef([data['U937_0.5_hrs'],data['U937_24_hrs']])
U937_comparison_2v4corr = U937_compar_2v4[0,1]


U937_compar_3v4=
np.corrcoef([data['U937_4_hrs'],data['U937_24_hrs']])
U937_comparison_3v4corr = U937_compar_3v4[0,1]



b= (('0v0.5', U937_comparison_1v2corr),('0v4', U937_comparison_1v3corr),'0v24', U937_comparison_1v4corr,'0.5v4', U937_comparison_2v3corr,'4v24', U937_comparison_2v4corr,'4vs24', U937_comparison_3v4corr))


print max(b)
"""
# PART (D)
np.sort(abs(data['Jurkat_0_hrs']-data['Jurkat_24_hrs call']))


### James' attempt at (D)
genelist = [['gene','value']]
for x in range(1,data_nodup.shape[0]):
    genelist = genelist + [[data['Gene Description'][x],abs(data['Jurkat_0_hrs'][x]-data['Jurkat_24_hrs call'][x])]]

genelist.sort(key = lambda x: x[::-1])
print genelist[1][0]
print genelist[2][0]
print genelist[3][0]
print genelist[4][0]
print genelist[5][0]
print genelist[6][0]
print genelist[7][0]
print genelist[8][0]
print genelist[9][0]
print genelist[10][0]




# PART (E)
#go through the list and check if any gene expression increase 2 fold at 24 hrs compared to 0 for all of the cell types. 
#if 24hr/0hr >=2 ---> add the gene to a list

answerE=range(0)
for gene in range(len(data)):
    if data.iat[gene,8]/data.iat[gene,2]>=2 and data.iat[gene,16]/data.iat[gene,10]>=2 and data.iat[gene,22]/data.iat[gene,18]>=2 and data.iat[gene,34]/data.iat[gene,28]>=2:
        answerE.append(data.iat[gene,0])
        
print answerE
    
# Part F:
#create a ansF.txt file

script, filename = argv


target = open(filename, 'w')
target.truncate()
for gene in range(len(data)):
	if data.iat[gene,2]/data.iat[gene,10]>=2.0 or data.iat[gene,2]/data.iat[gene,10] <=0.5:
	    target.write(data.iat[gene,1])
	    target.write('\n')
target.close()









