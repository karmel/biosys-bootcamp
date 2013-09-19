# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 18:13:57 2013

@author: Jeffrey
"""

from pandas import DataFrame,read_csv,Series
from pandas.io import parsers
import numpy

class ExpressionAnalyzer(object):
    '''
    Imports and analyzes gene expression data.
    Expects sample data as columns and genes as rows.
    '''
    
    def import_file(self, filename):
        '''Imports data from filename and returns a DataFrame object.'''
        return read_csv(filename)
    
    def num_genes(self, data):
        '''Prints out the number of unique 'Gene Description' values in data'''
        gene_set = set(data['Gene Description'])
        print len(gene_set)
    
    def correlated_data(self, data):
        '''Compares time points for each cell type and prints out the two most
           similar time points for each cell type''' 
        cell_types = [['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs'],
                      ['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs'],
                      ['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs'],
                      ['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']]
        
        for t in range(0,len(cell_types)):
            diffs = {}
            #Compares each time point with every other time point per cell type
            for i in range(0,len(cell_types[t])-1):
                for j in range(i+1,len(cell_types[t])):
                    #Returns the abs of the mean of the differences in values
                    diffs[abs(mean(data[cell_types[t][i]]-data[cell_types[t][j]]))] = cell_types[t][i].split("_")[1]+" "+cell_types[t][j].split("_")[1]
            #Uses a dict key to access the time point label of the result
            print cell_types[t][0].split("_")[0],diffs[min(diffs.keys())]
            
    def similar_cell(self, data):
        '''Aggregates data by cell type and prints out the data description 
           for a visual comparison'''
        data['HL'] = data['HL60_0_hrs']+data['HL60_0.5_hrs']+data['HL60_4_hrs']+data['HL60_24_hrs']
        data['U9'] = data['U937_0_hrs']+data['U937_0.5_hrs']+data['U937_4_hrs']+data['U937_24_hrs']
        data['NB'] = data['NB4_0_hrs']+data['NB4_5.5_hrs']+data['NB4_24_hrs']+data['NB4_48_hrs']+data['NB4_72_hrs']
        data['JU'] = data['Jurkat_0_hrs']+data['Jurkat_0.5_hrs']+data['Jurkat_4_hrs']+data['Jurkat_24_hrs']
        print 'HL60',data['HL'].describe(),'\n'
        print 'U937',data['U9'].describe(),'\n'
        print 'NB4',data['NB'].describe(),'\n'
        print 'Jurkat',data['JU'].describe(),'\n\n'
    
    def calibration_candidate(self, data):
        '''Finds the standard deviation of each cell type over time per gene
           and prints the ten genes with the lowest standard deviations'''
        std_list = []
        cell_types = [['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs'],
                      ['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs'],
                      ['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs'],
                      ['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']]
        
        for idx, row in data.iterrows():
            stdev = []
            for i in range(0,len(cell_types)):
                dat = []
                #Makes a list for each cell type
                for j in range(0,len(cell_types[i])):
                    dat.append(row[cell_types[i][j]])
                #Finds the standard deviations for each cell type per gene
                stdev.append(numpy.std(numpy.array([dat])))
            #Finds the means of the cell type standard deviations per gene
            std_list.append((data.ix[idx][1],mean(stdev)))
        #Returns the ten genes with the lowest standard deviations    
        print sorted(std_list, key=lambda x: x[1])[:10]

    def higher_expression_genes(self, data):
        '''Prints out a list of genes for which the 24-hr time point is
           two-fold higher than the 0-hr time point for all cell types'''
        hegenes = []
        cell_types = [['HL60_0_hrs','HL60_24_hrs'],
                      ['U937_0_hrs','U937_24_hrs'],
                      ['NB4_0_hrs','NB4_24_hrs'],
                      ['Jurkat_0_hrs','Jurkat_24_hrs']]
        
        for idx, row in data.iterrows():
            all_higher = True
            for i in range(0,len(cell_types)):
                    if row[cell_types[i][0]]*2 > row[cell_types[i][1]]:
                        all_higher = False
            if all_higher:
                hegenes.append(data.ix[idx][1])
        print '\n',hegenes
    
    def diff_reg_genes(self, data):
        '''Prints out a list of all genes differentially regulated in HL60
           cells as compared to U937 cells at 0 hrs'''
        drgenes = []
        for idx, row in data.iterrows():
            hl = row['HL60_0_hrs']
            u = row['U937_0_hrs']

            if hl >= u*2 or hl*2 <= u:
                drgenes.append(data.ix[idx][1])
        print '\n\n',drgenes
        
    
if __name__ == '__main__':
    ea = ExpressionAnalyzer()
    data = ea.import_file('data_set_HL60_U937_NB4_Jurkat.csv')
    ea.num_genes(data)
    ea.correlated_data(data)
    ea.similar_cell(data)
    ea.calibration_candidate(data)
    ea.higher_expression_genes(data)
    ea.diff_reg_genes(data)