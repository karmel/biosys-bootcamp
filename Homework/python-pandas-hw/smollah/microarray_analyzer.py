from pandas import DataFrame
from pandas.io import parsers
import pandas
from scipy.stats.stats import *

class ExpressionAnalyzer(object):
'''
Responsible for importing and analyzing gene expression data.
Expects data with named samples as columns and genes as rows.
'''
def import_file(self, filename):
''' Import data from file and return a DataFrame object. '''
    df = pandas.read_csv(filename)
    return df

def uniq_gene(df):
''' return unique genes. '''
    return len(set(df['Gene Description']))

def corr_matrix(df):
''' return pairwise correlation matrix. '''
    return df.corr()


if __name__ == '__main__':
ea = ExpressionAnalyzer()
ea.import_file(''/Users/smollah/Desktop/data_set_HL60_U937_NB4_Jurkat2.csv'')
