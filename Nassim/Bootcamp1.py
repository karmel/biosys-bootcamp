import pandas
import numpy as np
import math

#Note, this is incomplete, probably wrong, and certainly ugly.

class ExpressionAnalyzer(object):
    '''
    Responsible for importing and analyzing gene expression data.
    Expects data with named samples as columns and genes as rows.
    '''
    
    def __init__(self):
        pass
        
    def import_file(self, filename):
        '''
        Imports data from file and returns as DataFrame object
        '''
        
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()
        
        return df
 
 
 
        
    def DistinctGenes(self, filename):
        #Problem A
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()
        #The length of any given column should give the number of distinct genes
        print len(df['call'])



    
    def HighestCorrelation(self, filename):
        #Problem B
        
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()
        
        #Euclidean distance between two given values
        def eu_dist(name1, name2):
            return np.linalg.norm(name1 - name2)
            
        Time = [0, 0.5, 1, 24]
        HL60 = ['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs']
        test = float('Inf')
        
        #Check each time point combination, make sure they're unique
        for i in range(0,len(HL60)):
            for j in range(0,len(HL60)):
                if i < j:
                    compare = eu_dist(df[HL60[i]],df[HL60[j]])
        
                    #If the currently compared value is closer than the test value, it becomes the new "best" test value
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most highly correlated time points are "+str(Time[besti]) + " and " + str(Time[bestj])
        print "for HL60, at a euclidean distance value of " + str(test) + "." 

        #Repeat this block of code for each cell type
        U937 = ['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs']
        test = float('Inf')
        for i in range(0,len(U937)):
            for j in range(0,len(U937)):
                if i < j:
                    compare = eu_dist(df[U937[i]],df[U937[j]])
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most highly correlated time points are "+str(Time[besti]) + " and " + str(Time[bestj])
        print "for U937, at a euclidean distance value of " + str(test) + "."       

        Time2 = [0,5.5,24,48,72]
        NB4  = ['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']
        test = float('Inf')
        for i in range(0,len(NB4)):
            for j in range(0,len(NB4)):
                if i < j:
                    compare = eu_dist(df[NB4[i]],df[NB4[j]])
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most highly correlated time points are "+str(Time2[besti]) + " and " + str(Time2[bestj])
        print "for NB4, at a euclidean distance value of " + str(test) + "."    
        
        Jurkat=['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']
        test = float('Inf')
        for i in range(0,len(Jurkat)):
            for j in range(0,len(Jurkat)):
                if i < j:
                    compare = eu_dist(df[Jurkat[i]],df[Jurkat[j]])
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most highly correlated time points are "+str(Time[besti]) + " and " + str(Time[bestj])
        print "for Jurkat, at a euclidean distance value of " + str(test) + "."
        
        
    def SimilarCellTypes(self, filename):
        #Problem C
        
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()

        #Euclidean distance between two given values
        def eu_dist(name1, name2):
            return np.linalg.norm(name1 - name2)
        
        #Same concept as solution B, only instead of distinguishing by time points, we distinguish by different cells
        ZeroHour = ['HL60_0_hrs','U937_0_hrs','NB4_0_hrs','Jurkat_0_hrs']
        TwentyFourHour = ['HL60_24_hrs','U937_24_hrs','NB4_24_hrs','Jurkat_24_hrs']
        Names = ['HL60','U937','NB4','Jurkat']
        test = float('Inf')
        for i in range(0,len(ZeroHour)):
            for j in range(0,len(ZeroHour)):
                if i < j:
                    compare = eu_dist(df[ZeroHour[i]],df[ZeroHour[j]])
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most similar cell types are "+str(Names[besti]) + " and " + str(Names[bestj])
        print "at the 0th hour, at a euclidean distance value of " + str(test) + "."        
        
        test = float('Inf')
        for i in range(0,len(TwentyFourHour)):
            for j in range(0,len(TwentyFourHour)):
                if i < j:
                    compare = eu_dist(df[TwentyFourHour[i]],df[TwentyFourHour[j]])
                    if compare < test:
                        test = compare
                        besti = i
                        bestj = j
        print "The two most similar cell types are "+str(Names[besti]) + " and " + str(Names[bestj])
        print "at the 24th hour, at a euclidean distance value of " + str(test) + "."
        
        
    def GeneCandidates(self,filename):     
        #Problem D
        '''
        Find the ten genes with the smallest standard deviation between cell types, at t = 0 hr
        '''
        
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()
        
        Length = len(df['call']) - 1
        ZeroHour = ['HL60_0_hrs','U937_0_hrs','NB4_0_hrs','Jurkat_0_hrs']
        
        Temp = []
        StdDev = 0
        
        #These array values will be continuously replaced by smaller and smaller numbers
        MinValues = [10000000000,100000000001,100000000002,100000000003,100000000004,100000000005,100000000006,100000000007,\
                    100000000008,100000000009]
        
        #This array helps us correlate the smallest values 
        MinIdx = [0]*10
        
        for j in range(0,Length):
            for i in range(0,len(ZeroHour)):
                Temp.append(df[ZeroHour[i]].ix[j])
            StdDev = (np.std(Temp,axis=0))
            if max(MinValues) > StdDev:
                MinValues[MinValues.index(max(MinValues))] = StdDev
                MinIdx[MinValues.index(max(MinValues))] = j
        
        #print MinValues
        #print MinIdx
        #print min(MinValues)
        #print MinValues.index(min(MinValues))
        #print MinIdx[MinValues.index(min(MinValues))]
        print "Best conserved genes among cells:"
        for i in range(0,len(MinIdx)):
            print df['Gene Accession Number'].ix[MinIdx[i]]
        #Store = df['Gene Accession Number'].ix[MinIdx[MinValues.index(min(MinValues))]]
        #print "The best-conserved gene among cells is the " + str(Store) + " gene."
 
        
    def TwoFold(self,filename):     
        #Problem E
        
        '''This isn't done, and certainly not correct.'''
        
        data = open(filename,'r')
        df = pandas.read_csv(data,index_col=0,sep='\t')
        data.close()
        
        filtered = df[(2*df['HL60_0_hrs'] < df['HL60_24_hrs'])]
        filtered1 = filtered[(2*df['U937_0_hrs'] < df['U937_24_hrs'])]
        filtered2 = filtered1[(2*df['NB4_0_hrs'] < df['NB4_24_hrs'])]
        filtered3 = filtered2[(2*df['Jurkat_0_hrs'] < df['Jurkat_24_hrs'])]
        
        print filtered['Gene Accession Number']
        
        
        
    
                
        
        
            
        
        
        
            
            
        
    
