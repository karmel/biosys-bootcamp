
from __future__ import division
from pandas import DataFrame
from pandas.io import parsers
import re
import numpy as np

##
#Funtions#
##

#class ExpressionAnalyzer(object):#
#    '''Responsible for importing and analyzing gene expression data.'''#
#    '''Expects data with named samples as columns and genes as rows.'''#

#Function that imports a data file and return a DataFrame object.#
#For importing the data, it takes into account the rows that will be skipped and if the file contains a header#
def import_file(filename, skip_rows, header_pos):
    dataframe=parsers.read_csv(filename, skiprows=skip_rows, header=header_pos)
    return dataframe

#Function that counts the number of genes present in the file.#
#It takes into account that there are some gene ids of the form NN(...)NN_x, where NN(...)NN and NN(...)NN_x may be isoforms of the same gene and should be counted only once#
def genes_number(dataframe):
    num_info_col=1 #Number of the column that contains the gene ids#
    number=[]
    #num_info_col=1#
    for i in range(0,len(dataframe[dataframe.columns[num_info_col]])):
        number.append(re.sub(r'\_\w+', '', dataframe[dataframe.columns[num_info_col]][i]))
    
    return len(list(set(number)))

#Function that calculates for a given cell type (which has to be especified in the parameters), the correlation between its time points#
def time_point_corr(dataframe, num_info_col):
    dataframe_corr=dataframe[dataframe.columns[num_info_col]].corr()
    return dataframe_corr
#I will measure the similarity as the average of the correlation values that are obtained when comparing each pair of cell types throughout the diferent time points#


#Function that determines the gene's expression coefficient of variation for all the cell types during time zero, and based on that produces a list of the genes with the lowest coefficient of variation, those that should vary very little across samples#
def candidate_genes(dataframe, num_info_col, num_info_cols, candidate_num):
    coeff_var=[]
    coeff_var_abs=[]
    for i in range(0,len(dataframe[dataframe.columns[num_info_col]])):
        coeff_var.append(dataframe[dataframe.columns[num_info_cols]].ix[i].std()/dataframe[dataframe.columns[num_info_cols]].ix[i].mean())
        coeff_var_abs.append(abs(coeff_var[i]))
    
    gene_id=[]
    coeff_sorted=sorted(coeff_var_abs)[0:candidate_num]
    for i in range(0,len(coeff_sorted)):
        gene_id.append(coeff_var_abs.index(coeff_sorted[i]))
    
    genes_list=dataframe[dataframe.columns[num_info_col]].ix[gene_id]
    return genes_list

#Function that produces a list with the genes that show at least a two-fold higher expression when comparing two time points in all cell types#
def two_fold_higher(dataframe, num_info_col, num_cols_cond1, num_cols_cond2):
    two_fold_list=[]
    for i in range(0,len(dataframe[dataframe.columns[num_info_col]])):
        k=0
        for j in range(0, len(num_cols_cond1)):
            if ((dataframe[dataframe.columns[num_cols_cond2[j]]].ix[i]-dataframe[dataframe.columns[num_cols_cond1[j]]].ix[i])/dataframe[dataframe.columns[num_cols_cond1[j]]].ix[i])>=2:
                k=k+1
                if k==4: #This means that I will keep record of the gene only if the four cell types have a higher expression#
                    two_fold_list.append(dataframe[dataframe.columns[num_info_col]].ix[i])
            else:
                break
        
    
    return two_fold_list

#Function that produces an output file with the genes that show at least two-fold higher or lower expression, when comparing two cell types (especified as parameters)#
def two_fold_higher_lower(dataframe, num_info_col, col_cond1, col_cond2, output_file):
    two_fold_hl_list=[]
    for i in range(0,len(dataframe[dataframe.columns[num_info_col]])):
        if ((dataframe[dataframe.columns[col_cond2]].ix[i]-dataframe[dataframe.columns[col_cond1]].ix[i])/dataframe[dataframe.columns[col_cond1]].ix[i])>=2 or ((dataframe[dataframe.columns[col_cond2]].ix[i]-dataframe[dataframe.columns[col_cond1]].ix[i])/dataframe[dataframe.columns[col_cond1]].ix[i])<=-2:
            two_fold_hl_list.append(dataframe[dataframe.columns[num_info_col]].ix[i])
    
    DataFrame(two_fold_hl_list).to_csv(output_file, sep=',', index=False, header=False)
    return 'Your file was created with %s genes' % len(two_fold_hl_list)

##
#Execution of the functions#
##

#Load the DataFrame#
dfi=import_file('/Users/IFLM/Desktop/data_set_HL60_U937_NB4_Jurkat.csv', [1,2], 0)

#Calculate the number of genes#
gn=genes_number(dfi)
gn

#Calculate the time correlations for all the cell types#
HL60=time_point_corr(dfi, [2,4,6,8])
HL60
U937=time_point_corr(dfi, [10,12,14,16])
U937
NB4=time_point_corr(dfi, [18, 20, 22, 24, 26])
NB4
Jurkat=time_point_corr(dfi, [28, 30, 32, 34])
Jurkat

#Calculate which  cell types are the most similar#
HL60_U937=time_point_corr(dfi, [2,10,18,28]).ix[0,1], time_point_corr(dfi, [4,12,30]).ix[0,1], time_point_corr(dfi, [6,14,20,32]).ix[0,1], time_point_corr(dfi, [8,16,22,34]).ix[0,1]
HL60_NB4=time_point_corr(dfi, [2,10,18,28]).ix[0,2], time_point_corr(dfi, [6,14,20,32]).ix[0,2], time_point_corr(dfi, [8,16,22,34]).ix[0,2]
HL60_Jurkat=time_point_corr(dfi, [2,10,18,28]).ix[0,3], time_point_corr(dfi, [4,12,30]).ix[0,2], time_point_corr(dfi, [6,14,20,32]).ix[0,3], time_point_corr(dfi, [8,16,22,34]).ix[0,3]
U937_NB4=time_point_corr(dfi, [2,10,18,28]).ix[1,2], time_point_corr(dfi, [6,14,20,32]).ix[1,2], time_point_corr(dfi, [8,16,22,34]).ix[1,2]
U937_Jurkat=time_point_corr(dfi, [2,10,18,28]).ix[1,3], time_point_corr(dfi, [4,12,30]).ix[1,2], time_point_corr(dfi, [6,14,20,32]).ix[1,3], time_point_corr(dfi, [8,16,22,34]).ix[1,3]
NB4_Jurkat=time_point_corr(dfi, [2,10,18,28]).ix[2,3], time_point_corr(dfi, [6,14,20,32]).ix[2,3], time_point_corr(dfi, [8,16,22,34]).ix[2,3]
#I will measure the similarity as the average of the correlation values that are obtained when comparing each pair of cell types throughout the diferent time points#
np.array(HL60_U937).mean()
np.array(HL60_NB4).mean()
np.array(HL60_Jurkat).mean()
np.array(U937_NB4).mean()
np.array(U937_Jurkat).mean()
np.array(NB4_Jurkat).mean()

#Search for the genes that vary very little across samples#
candidate_genes(dfi, 1, [2,10,18,28], 10)

#Search for the genes that show two-fold higher expression at 24 hours versus 0 hours for all four cell types#
two_fold_higher(dfi, 1, [2,10,18,28], [8,16,22,34])

#Search for the genes that show at least two-fold higher or lower expression in HL60 cells compared to U937 cells at 0 hours; it produces an output file with the list of genes.#
two_fold_higher_lower(dfi, 1, 10, 2, 'output_file.txt')
   

