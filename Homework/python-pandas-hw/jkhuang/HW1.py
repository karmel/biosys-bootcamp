	#-------------------------------------------------------------------------------
# Name:        HW1
# Purpose:
#
# Author:      Justin
#
# Created:     16/09/2013
# Copyright:   (c) Justin 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import pandas as pd
import itertools
import numpy
pd.set_option('display.width', 10000) #For ipython notebook

def QuestionA(expressiondata):
    print 'A. How many distinct genes are represented in the data set?'
    print len(set(expressiondata.index))

def QuestionB(expressiondata):
    print 'B. Which two time points are the most highly correlated for each cell type?'
    HL60 = ['HL60_0_hrs', 'HL60_0.5_hrs', 'HL60_4_hrs', 'HL60_24_hrs']
    U937 = ['U937_0_hrs', 'U937_0.5_hrs', 'U937_4_hrs', 'U937_24_hrs']
    NB4 = ['NB4_0_hrs', 'NB4_5.5_hrs', 'NB4_24_hrs', 'NB4_48_hrs', 'NB4_72_hrs']
    Jurkat = ['Jurkat_0_hrs', 'Jurkat_0.5_hrs', 'Jurkat_4_hrs', 'Jurkat_24_hrs']

    cellLines = [HL60, U937, NB4, Jurkat]
    for cellLine in cellLines:
        correlations=[]
        for timepoints in list(itertools.combinations(cellLine, 2)):
            correlations.append((timepoints[0], timepoints[1], abs(mean(expressiondata[timepoints[1]]-expressiondata[timepoints[0]]))))
        print sorted(correlations, key=lambda x: x[2])[0]

def QuestionC(expressiondata):
    print 'C. Which two cell types are the most similar?'
    #Comparing differentiated cell lines at 24 hours
    cellLines = ['HL60_24_hrs', 'U937_24_hrs', 'NB4_24_hrs', 'Jurkat_24_hrs']

    comparisons = []
    for combo in list(itertools.combinations(cellLines, 2)):
        comparisons.append((combo[0], combo[1], abs(mean(expressiondata[combo[1]]-expressiondata[combo[0]]))))
    for comparison in sorted(comparisons, key=lambda x: x[2]):
        print comparison

def QuestionD(expressiondata):
    print 'D. It is often useful to know which genes change very little across samples for the sake of normalization or calibration. Based on this data set, what are ten good candidates for genes to use to calibrate machinery or analyses across all these samples?'
    changes = []
    for i in range(len(expressiondata.index)):
        HL60expression = [expressiondata.ix[i]['HL60_0_hrs'], expressiondata.ix[i]['HL60_0.5_hrs'], expressiondata.ix[i]['HL60_4_hrs'], expressiondata.ix[i]['HL60_24_hrs']]
        U937expression = [expressiondata.ix[i]['U937_0_hrs'], expressiondata.ix[i]['U937_0.5_hrs'], expressiondata.ix[i]['U937_4_hrs'], expressiondata.ix[i]['U937_24_hrs']]
        NB4expression = [expressiondata.ix[i]['NB4_0_hrs'], expressiondata.ix[i]['NB4_5.5_hrs'], expressiondata.ix[i]['NB4_24_hrs'], expressiondata.ix[i]['NB4_48_hrs'], expressiondata.ix[i]['NB4_72_hrs']]
        Jurkatexpression = [expressiondata.ix[i]['Jurkat_0_hrs'], expressiondata.ix[i]['Jurkat_0.5_hrs'], expressiondata.ix[i]['Jurkat_4_hrs'], expressiondata.ix[i]['Jurkat_24_hrs']]

        HL60std = abs(numpy.std(numpy.array(HL60expression)))
        U937std = abs(numpy.std(numpy.array(U937expression)))
        NB4std = abs(numpy.std(numpy.array(NB4expression)))
        Jurkatstd = abs(numpy.std(numpy.array(Jurkatexpression)))
        changes.append((i, abs(mean((HL60std, U937std, NB4std, Jurkatstd)))))

    for item in list(sorted(changes, key=lambda x: x[1])[:10]):
        print expressiondata.index[item[0]], '---Avg Standard Deviations: ',item[1]

def QuestionE(expressiondata, AccessionNumber=False):
    print 'E. Do any genes show two-fold higher expression at 24 hours versus 0 hours for all four cell types? If so, which ones?'
    highexpressiongenes = []
    for i in range(len(expressiondata.index)):
        cellexpressions=[0,0,0,0]
        diff1 = float(expressiondata.ix[i]['HL60_24_hrs']) > 2 * float(expressiondata.ix[i]['HL60_0_hrs'])
        if diff1:
            cellexpressions[0]=1
        diff2 = float(expressiondata.ix[i]['U937_24_hrs']) > 2 * float(expressiondata.ix[i]['U937_0_hrs'])
        if diff2:
            cellexpressions[1]=1
        diff3 = float(expressiondata.ix[i]['NB4_24_hrs']) > 2 * float(expressiondata.ix[i]['NB4_0_hrs'])
        if diff3:
            cellexpressions[2]=1
        diff4 = float(expressiondata.ix[i]['Jurkat_24_hrs']) > 2 * float(expressiondata.ix[i]['Jurkat_0_hrs'])
        if diff4:
            cellexpressions[3]=1
        if sum(cellexpressions)==4:
            if AccessionNumber:
                highexpressiongenes.append(expressiondata.ix[i]['Gene Accession Number'])
            else:
                highexpressiongenes.append(expressiondata.index[i])

    print len(set(highexpressiongenes)), 'genes'
    for gene in set(highexpressiongenes): print gene

def QuestionF(expressiondata, AccessionNumber=False):
    print 'F. Which genes are differentially regulated (at least two-fold higher or lower) in HL60 cells as compared to U937 cells at 0 hours?'
    diffregulatedgenes = []
    for i in range(len(expressiondata.index)):
        diff1 = float(expressiondata.ix[i]['HL60_0_hrs']) > 2 * float(expressiondata.ix[i]['U937_0_hrs']) #HL60 two-fold higher expression than U937 @0 hours
        diff2 = float(expressiondata.ix[i]['U937_0_hrs']) > 2 * float(expressiondata.ix[i]['HL60_0_hrs']) #U937 two-fold higher expression than HL60 @0 hours
        if diff1:
            if AccessionNumber:
                diffregulatedgenes.append(expressiondata.ix[i]['Gene Accession Number'])
            else:
                diffregulatedgenes.append(expressiondata.index[i])
        elif diff2:
            if AccessionNumber:
                diffregulatedgenes.append(expressiondata.ix[i]['Gene Accession Number'])
            else:
                diffregulatedgenes.append(expressiondata.index[i])

    print len(set(diffregulatedgenes)), 'genes'
    for gene in set(diffregulatedgenes): print gene

def main():
    #New Data posted by Karmel to Github repository
    data = pd.read_csv('C:\Users\Justin\Documents\GitHub\\biosys-bootcamp\Homework\python-pandas-hw\data_set_HL60_U937_NB4_Jurkat.txt', sep='\t')
    cols = list(data.columns[1:])
    cols[-1]='Jurkat_24_hrs'
    cols.append('call.16')
    data.columns=cols #Fixing last column alignment

    QuestionA(data)
    QuestionB(data)
    QuestionC(data)
    QuestionD(data)
    QuestionE(data)
    QuestionF(data)

if __name__ == '__main__':
    main()
