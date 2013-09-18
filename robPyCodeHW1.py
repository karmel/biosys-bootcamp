# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

cd ~/Documents/Projects/PythonHomework/

# <codecell>

import numpy as np
from pandas import *
from scipy import stats
from __future__ import division

class ExpressionTableAnalyzer(object):
    ''' Allows the analyzer of dataFrame objects that contain microarray expression
    data '''
    
    
    def import_file(self, filename):
        ''' Imports a tab separated file and parses it into a pandas dataFrame. '''
        open_file = open(filename,'r')
        self.dataFrame = read_csv(open_file, sep='\t')

    def cfVariation(self, vector):
        '''Calculates the absolute value of the coefficient of variation from a
           vector or series input.'''
        return abs(stats.variation(vector, axis=0))

    
    def twoFoldHigherCalc(self, dataFrame, col_name1, col_name2):
        ''' Returns the row index of genes with 2 fold expression difference.In it's
current implentation there are false negatives which are not reported for
genes with expression which goes from negative to positive values. '''
        numRows = 7229
        iterRow = dataFrame[[col_name1, col_name2]]
        index_container = []
        for i in range(numRows):
            rowIteration = iterRow.ix[i]
            if ((rowIteration[0]<0 and rowIteration[0] > -1) and rowIteration[1] > 0):
                index_container.insert(0,-1)
            if ((rowIteration[0]/rowIteration[1])>2 and rowIteration[0]<0) or (((rowIteration[0]/rowIteration[1])<0.5 and ((rowIteration[0]/rowIteration[1])>0)) and (rowIteration[0]>0)):
                index_container.insert(0,i)
        return index_container
    
    def twoFoldCalc(self, dataFrame, col_name1, col_name2):
        ''' Returns the row index of genes with 2 fold expression difference.In it's
current implentation there are false negatives which are not reported for
genes with expression which goes from negative to positive values. '''
        numRows = 7229
        iterRow = dataFrame[[col_name1, col_name2]]
        index_container = []
        for i in range(numRows):
            rowIteration = iterRow.ix[i]
            if ((rowIteration[0]/rowIteration[1] > 2) or (rowIteration[0]/rowIteration[1] < 0.5 and rowIteration[0]/rowIteration[1]>0)):
                index_container.insert(0,i)
        return index_container

    def scoreExpressionVectors(self, dataFrame, columnName1, columnName2):
        ''' Returns the norm of the difference of two expression vector columns.'''
        column1 = dataFrame[columnName1]
        column2 = dataFrame[columnName2]
        return np.linalg.norm(column1-column2)


    def highScoreAllPermutations(self, dataFrame, listNames):
        ''' Calculates the similarity score between timepoints for every possible
combination given some list input of column header names for each timepoint to be compared. '''
        highScore = float('Inf')
        while len(listNames)>1:
            for i in range(len(listNames)-1):
                currentScore = self.scoreExpressionVectors(dataFrame, listNames[0], listNames[i+1])
            
                if currentScore<highScore:
                    highScorers = [listNames[0], listNames[i+1], currentScore]
                    highScore = currentScore
                
            listNames = listNames[1:]
        return highScorers
    
    def scoreCells(self, dataFrame, listCellTypes):
        lowScore = float('Inf')
        blue = 2
        while len(listCellTypes)>1:
            for i in range(blue):
                currentScore = scoreRow(dataFrame, listCellTypes[0], listCellTypes[i+1])
                if currentScore <= lowScore:
                    winnerRecord = [listCellTypes[0], listCellTypes[i+1], currentScore]
                    lowScore = currentScore
            listCellTypes = listCellTypes[1:]
            blue = blue-1
        return winnerRecord
    
    def scoreRow(self, dataFrame, colSubset1, colSubset2):
        rowSubset1 = dataFrame[colSubset1]
        rowSubset2 = dataFrame[colSubset2]
    
        rowScores1 = []
        rowScores2 = []
    
        for i in range(7229):
            geneCompare1 = rowSubset1.ix[i]
            geneCompare2 = rowSubset2.ix[i]
            
            rowDiff = range(4)
            
            for j in range(4):
                rowDiff[j] = geneCompare1[j] - geneCompare2[j]
            
            rowScores1.insert(0,np.linalg.norm(rowDiff))
        return sum(rowScores1)
    

if __name__ == '__main__':
    df = ExpressionTableAnalyzer()
    df.import_file('/Users/rkf/Documents/Projects/PythonHomework/data_set.txt')

    HL_cells = ['HL60_0_hrs', 'HL60_0.5_hrs', 'HL60_4_hrs', 'HL60_24_hrs']

    U937_cells = ['U937_0_hrs', 'U937_0.5_hrs', 'U937_4_hrs', 'U937_24_hrs']

    NB4_cells = ['NB4_0_hrs', 'NB4_5.5_hrs', 'NB4_24_hrs', 'NB4_48_hrs', 'NB4_72_hrs']

    Jurkat_cells = ['Jurkat_0_hrs', 'Jurkat_0.5_hrs', 'Jurkat_4_hrs', 'Jurkat_24_hrs']

    q_a = [HL_cells, U937_cells, NB4_cells, Jurkat_cells]

    tp0 = ['HL60_0_hrs', 'U937_0_hrs', 'Jurkat_0_hrs']

    tp5 = ['HL60_0.5_hrs', 'U937_0.5_hrs', 'Jurkat_0.5_hrs']

    tp4 = ['HL60_4_hrs', 'U937_4_hrs', 'Jurkat_4_hrs']

    tp24 = ['HL60_24_hrs', 'U937_24_hrs', 'NB4_24_hrs', 'Jurkat_24_hrs']

    print 'A: The number of unique genes as calculated by unique gene descriptions is'
    print(len(set(df.dataFrame['Gene Description'])))

    print 'B: The most similar correlated time points for each cell are '
    
    for i in q_a:
        print(df.highScoreAllPermutations(df.dataFrame,i))

    print 'C: Which two cell types are the most similar?'
    
    listOfLists = [HL_cells, Jurkat_cells, U937_cells]
    print df.scoreCells(df.dataFrame, listOfLists)

    tpPerCellType=df.dataFrame[HL_cells]
    variationHL = (tpPerCellType.apply(cfVariation, axis=1)).order()

    print 'D: 10 least variant genes by coefficient variation.'
    print df.dataFrame['Gene Accession Number'].ix[variationHL.index[0:9]]

    print 'E: Do any genes show two fold higher expression in all cell types?\nThe set below contains the row index of these genes.'
    winnerSetContainer = 0
    for i in q_a:
        counter =0
        tempVar = set(df.twoFoldHigherCalc(df.dataFrame, i[0], i[3]))
        if counter==0:
            winnerSetContainer=tempVar
        else:
            winnerSetContainer = winnerSetContainer.intersection(tempVar)
        counter +=1
    print winnerSetContainer
    
    print 'F: Which HL60 genes and U937 genes at zero hours are differentially expressed?'
    print df.dataFrame['Gene Accession Number'].ix[df.twoFoldCalc(df.dataFrame, HL_cells[0], U937_cells[0])]
    
    print 'G: Yo ho ho and a bottle of rum!'
    
    
        

# <codecell>


