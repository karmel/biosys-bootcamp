#!/home/jdavistu/epd-7.3-2-rh3-x86_64/bin/python

# Michael Kramer, Bioinformatics Bootcamp, 9/19/12

from __future__ import division
from pandas import DataFrame
from pandas.io import parsers
from scipy import stats

class ExpressionAnalyzer(object):
    '''
    Responsible for importing and analyzing gene expression data.  Expects data with named samples as columns and genes as rows
    '''
    
    def import_file(self, filename, separator):
        ''' Imports file '''
        self.data = parsers.read_csv(filename, sep=separator)

    def numUniqueGenes(self, geneListColName):
        ''' Returns the number of unique genes from the data frame when given a column name with unique gene identifiers'''
        return len(set(self.data[geneListColName]))

    def getMostCorrelatedTimePointsByCellType(self,cell_types):
        ''' Takes a list of cell types and returns a list of lists, where each element of the list has the form [cell_type, [time_point1 time_point2]] where time_point1 and time_point2 are the most correlated time points of that cell type by Pearson correlation'''
        # Descriptive comments! Yay!
        retList = []
        for cell_type in cell_types:
            max_cor = -2
            cell_type_expression_cols = [i for i in self.data.columns if i.find(cell_type) != -1]
            # List comprehension! Right on!
            for i in range(0,(len(cell_type_expression_cols)-1)):
                for j in range(i+1,len(cell_type_expression_cols)):
                    cor = stats.pearsonr(ea.data[cell_type_expression_cols[i]],ea.data[cell_type_expression_cols[j]])[0]
                    if cor > max_cor:
                        max_cor_pair = [cell_type_expression_cols[i], cell_type_expression_cols[j]]
                        max_cor = cor
            retList.append([cell_type,max_cor_pair])
        return retList

    def getMostSimilarCellTypes(self,cell_types,time_points):
        '''Takes a list of cell types and a list of time points to compare 
        and returns a list of two members which are the most similar cell types based on average Pearson correlation 
        accross the different time points'''
        # A general rule of style in Python is that lines should be less than 80 char long.
        # That's a little stringent and often not followed for the sake of clarity,
        # but if you're the kind of person who writes super long lines of code, 
        # it is worthwhile to practice breaking them up. Makes the code easier to read for us reviewers.
        max_sim = -2
        for i in range(0,(len(cell_types)-1)):
            for j in range(i+1,len(cell_types)):
                similarity = 0;
                num_same_time_points = 0;
                for time_point in time_points:
                    cell_type_1_tp = [k for k in ea.data.columns if ((k.find(cell_types[i]) != -1) and (k.find('_'+time_point+'_') != -1))]
                    cell_type_2_tp = [l for l in ea.data.columns if ((l.find(cell_types[j]) != -1) and (l.find('_'+time_point+'_') != -1))]
                    if ((len(cell_type_1_tp) == 1) and (len(cell_type_2_tp) == 1)):
                        cor = stats.pearsonr(ea.data[cell_type_1_tp[0]],ea.data[cell_type_2_tp[0]])[0]
                        similarity += cor
                        num_same_time_points += 1
                similarity /= num_same_time_points
                if similarity > max_sim:
                    max_sim = similarity
                    max_sim_pair = [cell_types[i],cell_types[j]]
        return max_sim_pair

    def getLeastVariantGenes(self,geneColumnName,numToGet):
        ''' Takes a column name of the column which contains gene identifiers and a number of genes to get and returns a list of the 10 least varying genes accross all conditions'''
        # Probably more looping than is necessary-- can be done with slices and no appending.
        geneList = []
        for i in self.data.var(axis=1).order()[0:numToGet-1].index:
            geneList.append(self.data[geneColumnName][i])
        return geneList

    def getGenesHigherInAllCellTypesAtTimePoint(self,initial_time,second_time,fold_change,cell_types,geneColumnName):
        ''' Takes an initial time, second time point, fold change, list of cell types to compare and a column name containing gene identifiers and returns a list of genes elevated in all cell_types at the second time point compared to the first time point by the given fold change'''
        geneList = []
        for gene_index in range(0,len(self.data.index)):
            add_gene = True
            for cell_type in cell_types:
                if self.data[cell_type+"_"+ second_time][gene_index] < fold_change*ea.data[cell_type+"_" + initial_time][gene_index]:
                    add_gene = False
                    break
            if (add_gene == True):
                geneList.append(ea.data[geneColumnName][gene_index])
        return geneList

    def getGenesDE(self,col1,col2,fold_change,geneColumnName):
        '''Takes two column names, a fold change and a column name containing gene identifiers and returns a list of genes differentially expressed between the two columns by at least the specified fold change'''
        # Love that you used a fold_change var. 
        geneList = []
        for gene_index in range(0,len(self.data.index)):
            if not (1/fold_change*ea.data[col1][gene_index] < ea.data[col2][gene_index] < fold_change*ea.data[col1][gene_index]):
                geneList.append(ea.data[geneColumnName][gene_index])
        return geneList
                
if __name__ == '__main__':
    
    # YAYYYY for following the template.
    print "Michael Kramer, Bioinformatics Bootcampe Homework 1, 9/19/2012\n"
    
    ea = ExpressionAnalyzer()
    ea.import_file('/home/mkramer/bootcamp/data_set_HL60_U937_NB4_Jurkat.txt',"\t")
    print "A) There are", ea.numUniqueGenes("Gene Accession Number"), "distinct genes"

    print "B) Most correlated time points by cell type: "
    cell_types = ['HL60', 'U937', 'NB4', 'Jurkat']
    max_cor_pairs = ea.getMostCorrelatedTimePointsByCellType(cell_types)
    for i in max_cor_pairs:
        print i[0], i[1]

    print "C) Most similar cell types: "
    time_points = ['0', '0.5', '4', '24']
    print ea.getMostSimilarCellTypes(cell_types,time_points)

    print "D) 10 genes with least variance: "
    ten_least_variant_genes = ea.getLeastVariantGenes("Gene Accession Number", 10)
    for i in  ten_least_variant_genes:
        print i

    # Okay, I tend towards really verbose var names too, but this might be a little excessive :)
    # See above about shorter lines so that reviewers don't have to scroll...
    genes_2_fold_higher_in_all_4_at_24_compared_to_0 = ea.getGenesHigherInAllCellTypesAtTimePoint('0_hrs','24_hrs',2,cell_types,"Gene Accession Number")
    print "E) Genes 2 fold higher at 24 hours compared to 0 hours in all 4 cell lines: "
    for gene in genes_2_fold_higher_in_all_4_at_24_compared_to_0:
        print gene
        

    genes_de_in_HL60_vs_U937 = ea.getGenesDE("U937_0_hrs","HL60_0_hrs",2,"Gene Accession Number")
    print "F) Genes differentially expressed in HL60 vs U937 cells at 0 hours: "
    for gene in genes_de_in_HL60_vs_U937:
        print gene
