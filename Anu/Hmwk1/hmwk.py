import numpy as np
import pandas
from pandas.io import parsers
from pandas import Series, DataFrame
import scipy as sp
import numpy as np
import re # was using this, but not anymore

# Importing File 
#reads in tab delimited data from a local file to turn into dataframe
#local file 'C:\Python27\data.txt'
#readin( 'C:\Python27\data.txt')
def readin (filename):
    data= pandas.io.parsers.read_table(filename,sep='\t*',skiprows=0 ,header= 0,na_values=' ')
    data= DataFrame(data)
    print data.columns
    return(data)

##Part A 
#Gets rid of redundant genes to tell you how many genes are in the array
#According to the paper there were 6416 genes in this dataset...but I keep getting 6570 distinct genes
#Run redundant(readin('C:\Python27\data.txt'),'Gene Description')
def redundant (df,item):
    a= df[item]
    a= a.fillna(0)
    print a 
    b=0 # counter for redundant gene descriptions
    c=0 # counter for nan
    for i in range(0,len(a)-1):
        if a[i] == 0 :
            c=c+1
        if a[i] == a[i+1]:
            b=b+1
    numgenes=((len(a))-b-c) # total genes-redundant gene descreptions- emptys 
    return (numgenes)
   


#Part B 
#Looks at the correlations between values for timepoints within a cell to determine which two time points are most similar for each of the four cell types
#run corwithincell(readin('C:\Python27\data.txt')

def corwithincell(df):
#this is to drop all the rows wtih non numeric values  
    df=df.drop('call', axis=1)
    df=df.drop('call.1',axis=1)
    df=df.drop('call.2',axis=1)
    df=df.drop('call.3',axis=1)
    df=df.drop('call.4',axis=1)
    df=df.drop('call.5',axis=1)
    df=df.drop('call.6',axis=1)
    df=df.drop('call.7',axis=1)
    df=df.drop('call.8',axis=1)
    df=df.drop('call.9',axis=1)
    df=df.drop('call.10',axis=1)
    df=df.drop('call.11',axis=1)
    df=df.drop('call.12',axis=1)
    df=df.drop('call.13',axis=1)
    df=df.drop('call.14',axis=1)
    df=df.drop('call.15',axis=1)
    df=df.drop('call.16',axis=1)

#HL60 cells at 4 timepoints  (should have appended all of these together, didn't have time to)
    cellone=df[df.columns[2]]
    cellonea=df[df.columns[3]]
    celloneb=df[df.columns[4]]
    cellonec=df[df.columns[5]]

#U937 cells at 4 timepoints
    celltwo=df[df.columns[6]]
    celltwoa=df[df.columns[7]]
    celltwob=df[df.columns[8]]
    celltwoc=df[df.columns[9]]

#NB4 cells at 5 timepoints 
    cellthree=df[df.columns[10]]
    cellthreea=df[df.columns[11]]
    cellthreeb=df[df.columns[12]]
    cellthreec=df[df.columns[13]]
    cellthreed=df[df.columns[14]]

#Jurkat cells at 4 timepoints  
    cellfour=df[df.columns[15]]
    cellfoura=df[df.columns[16]]
    cellfourb=df[df.columns[17]]
    cellfourc=df[df.columns[18]]

#Correlations for HL60 cells  
    a= cellone.corr(cellonea)
    b= cellone.corr(celloneb)
    c= cellone.corr(cellonec)
    d= cellonea.corr(celloneb)
    e= cellonea.corr(cellonec)
    f= celloneb.corr(cellonec)
    l=[a,b,c,d,e,f]
    val=l.index(max(l))
    tps = { '0and30min' :a, '0and4hrs':b, '0and24hrs':c,'30minand4hrs':d, '30 minand24hrs':e,'4hrsand24hrs':f }
    keys=tps.keys()
    print 'for HlL60 cells, the most correlated time points are:', print keys[val]

#Correlations for U937 cells
    a= celltwo.corr(celltwoa)
    b= celltwo.corr(celltwob)
    c= celltwo.corr(celltwoc)
    d= celltwoa.corr(celltwob)
    e= celltwoa.corr(celltwoc)
    f= celltwob.corr(celltwoc)
    l=[a,b,c,d,e,f]
    val=l.index(max(l))
    tps = { '0and30min' :a, '0and4hrs':b, '0and24hrs':c,'30minand4hrs':d, '30 minand24hrs':e,'4hrsand24hrs':f }
    keys=tps.keys()
    print 'for U937 cells, the most correlated time points are:',print keys[val]

# Correlations for NB4 cells
    a= cellthree.corr(cellthreea)
    b= cellthree.corr(cellthreeb)
    c= cellthree.corr(cellthreec)
    d= cellthree.corr(cellthreed)
    e= cellthreea.corr(cellthreeb)
    f= cellthreea.corr(cellthreec)
    g= cellthreea.corr(cellthreed)
    h= cellthreeb.corr(cellthreec)
    i= cellthreeb.corr(cellthreed)
    j= cellthreec.corr(cellthreed)
    l=[a,b,c,d,e,f,g,h,i,j]
    val=l.index(max(l))
    tps = { '0and5.5hrs' :a, '0and24hrs':b, '0and48hrs':c,'0and72hrs':d, '5.5hrsand24hrs':e,'5.5hrsand48hrs':f ,'5.5hrsand72hrs':g,'24hrsand48hrs':h,'24hrsand72hrs':i,'48hrsand72hrs':j}
    keys=tps.keys()
    print 'for NB4 cells, the most correlated time points are:',print keys[val]

# Correlations for Jurkat Cells
    a= cellfour.corr(cellfoura)
    b= cellfour.corr(cellfourb)
    c= cellfour.corr(cellfourc)
    d= cellfour.corr(cellfourb)
    e= cellfour.corr(cellfourc)
    f= cellfour.corr(cellfourc)
    l=[a,b,c,d,e,f]
    val=l.index(max(l))
    tps = { '0and30min' :a, '0and4hrs':b, '0and24hrs':c,'30minand4hrs':d, '30 minand24hrs':e,'4hrsand24hrs':f }
    keys=tps.keys()
    print 'for Jurkat cells, the most correlated time points are:'
    print keys[val]


##Part C
#Finds which two celltypes are the most similar by inputting tpforx(number of timepoints of data)
#and the number of timepoints taken at less than or equal to 24 hours for 4 celltypes over a period of 24 hours 
#Computed median value for each cell type across timepoints, and then found smallest difference amongst median values
#Dictionary names are preset to the four cell types in this data set
#run mostsimilar(readin('C:\Python27\data.txt'),4,0,4,0,5,3,4,0)

def mostsimilar(df,tpfor1,numbelow24,tpfor2,num1below24,tpfor3,num2below24,tpfor4,num3below24):
    medians=df.median(axis=0)
    medians= [np.median(medians[0:numbelow24-1]),np.median(medians[tpfor1:tpfor1+num1below24-1]),np.median(medians[tpfor2+tpfor1:tpfor1+tpfor2+num2below24-1]),np.median(medians[tpfor1+tpfor2+tpfor3:tpfor1+tpfor2+tpfor3+num3below24-1])]
    dif= [abs(medians[0]-medians[1]), abs(medians[0]-medians[2]), abs(medians[0]-medians[3]),abs(medians[1]-medians[2]),abs(medians[1]-medians[3]),abs(medians[2]-medians[3])]
    celltyps= {'HL60andU937':abs(medians[0]-medians[1]), 'HL60andNB4' :abs(medians[0]-medians[2]), 'HL60andJurkat':abs(medians[0]-medians[3]),'U937andNB4':abs(medians[1]-medians[2]),'U937andJurkat':abs(medians[1]-medians[3]), 'NB4andJurkat':abs(medians[2]-medians[3])}
    val= dif.index(min(dif))
    keys=celltyps.keys()
    print 'Cell types',keys[val],'are the most similar'


##Part D
#Which genes change very little across samples
#finds the top 10 genes with the smallest standard deviation over time/cell type
#controls(readin('C:\Python27\data.txt'),'Gene Description')
def controls(df,genedesc_item):
    sd=df.std(axis=1)
    sd=DataFrame(sd)
    newdata= df.join(sd)
    newdata= newdata.sort_index(axis=0,by= 0L,ascending= True)
    candidates= newdata[genedesc_item][2:12]
    return (candidates)



##Part E
#Function to show which genes are two fold upregulated at 24 hours 
#upreg(readin('C:\Python27\data.txt'),'Gene Description')
def upreg(df,genedesc_item):
#data= pandas.io.parsers.read_table('C:\Python27\data.txt',sep='\t*',skiprows=0 ,header= 0,na_values=' ')

#this is to drop all the rows wtih non numeric values  
    df=df.drop('call', axis=1)
    df=df.drop('call.1',axis=1)
    df=df.drop('call.2',axis=1)
    df=df.drop('call.3',axis=1)
    df=df.drop('call.4',axis=1)
    df=df.drop('call.5',axis=1)
    df=df.drop('call.6',axis=1)
    df=df.drop('call.7',axis=1)
    df=df.drop('call.8',axis=1)
    df=df.drop('call.9',axis=1)
    df=df.drop('call.10',axis=1)
    df=df.drop('call.11',axis=1)
    df=df.drop('call.12',axis=1)
    df=df.drop('call.13',axis=1)
    df=df.drop('call.14',axis=1)
    df=df.drop('call.15',axis=1)
    df=df.drop('call.16',axis=1)

#HL60 cells at 4 timepoints  (should have appended all of these together, didn't have time to)
    cellone=df[df.columns[2]]
    cellonea=df[df.columns[3]]
    celloneb=df[df.columns[4]]
    cellonec=df[df.columns[5]]

#U937 cells at 4 timepoints
    celltwo=df[df.columns[6]]
    celltwoa=df[df.columns[7]]
    celltwob=df[df.columns[8]]
    celltwoc=df[df.columns[9]]

#NB4 cells at 5 timepoints 
    cellthree=df[df.columns[10]]
    cellthreea=df[df.columns[11]]
    cellthreeb=df[df.columns[12]]
    cellthreec=df[df.columns[13]]
    cellthreed=df[df.columns[14]]

#Jurkat cells at 4 timepoints  
    cellfour=df[df.columns[15]]
    cellfoura=df[df.columns[16]]
    cellfourb=df[df.columns[17]]
    cellfourc=df[df.columns[18]]

    length=len (df[df.columns[1]])
    HL60=[]
    U937=[]
    NB4=[]
    Jurkat=[]
                 
    for i in range(0,length-1): #cheking for 2 fold upregulation
        if cellonec[i]==cellone[i]*2:
            HL60.append(df[genedesc_item][i])
        if celltwoc[i]==celltwo[i]*2:
            U937.append(df[genedesc_item][i])
        if cellthreec[i]==cellthree[i]*2:
            NB4.append(df[genedesc_item][i])
        if cellfourc[i]==cellfour[i]*2:
            Jurkat.append(df[genedesc_item][i])

#printing all the genes that are two fold upregulated
    print'In HL60 cells, these genes are upregulated 2fold at 24 hours', HL60
    print'In U937 cells, these genes are upregulated 2fold at 24 hours', U937
    print'In NB4 cells, these genes are upregulated 2fold at 24 hours', NB4
    print'In Jurkat cells, these genes are upregulated 2fold at 24 hours', Jurkat

#Writing Genes to CSV in case needed in the future (2 fold upregulated at 24 hours)
    HL60=DataFrame(HL60)
    HL60.to_csv('C:\Python27\HL60.csv',sep=',')
    U937=DataFrame(U937)
    U937.to_csv('C:\Python27\U937.csv',sep=',')
    NB4=DataFrame(NB4)
    NB4.to_csv('C:\Python27\NB4.csv',sep=',')
    Jurkat=DataFrame(Jurkat)
    Jurkat.to_csv('C:\Python27\Jurkat.csv',sep=',')

##Part F
# Printing genes two fold down regulated at t=0 for HL60 cs U937
# Saves list to csv to be uploaded into David
# run upreg(readin('C:\Python27\data.txt'),'Gene Description','Gene Accession Number')
def downreg(df,genedesc_item,geneid_item):
    cellone=df[df.columns[2]]
    celltwo=df[df.columns[6]]
    downreg=[]
    accession=[]
    for i in range(0,length-1):
        if (cellone[i]>=celltwo[i]*2):
            downreg.append(df[genedesc_item][i])
            accession.append(df[geneid_item][i])
        if (2*cellone[i]<= celltwo[i]):
            downreg.append(df[genedesc_item][i])
            accession.append(df[geneid_item][i])
        
    downreg=DataFrame(downreg)
    downreg.to_csv('C:\Python27\downreg.csv',sep=',')
    accession=DataFrame(accession)
    accession.to_csv('C:\Python27\list.csv',sep=',')
