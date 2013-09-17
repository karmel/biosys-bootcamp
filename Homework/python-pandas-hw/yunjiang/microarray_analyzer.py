#usage ipython microarray_analyzer.py >answers.txt
import pandas
from pandas import DataFrame
from pandas.io import parsers

class ExpressionAnalyzer(object):
	'''
	Responsible for importing and analyzing gene expression data.
	Expects data with named samples as columns and genes as rows.
	'''

	def _init_(self):
		self.data = {};

	def import_file(self, filename):
		''' Import data from file and return a DataFrame object. '''
		self.data = pandas.read_csv(filename,sep='\t',index_col=False)

	def count_gene(self):
		'''Count gene number'''
		print "The gene number is:"
		gene_set = set (self.data["Gene Accession Number"])
		print len(gene_set) 	
		print '\n'

	def time_corr(self):
		'''Compute correlation coefficient between time points for each cell type'''
		print 'The corelation matrix is:\n'	
		print self.data[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs']].corr()
		print '\n'
		print self.data[['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs']].corr()
		print '\n'
		print self.data[['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs call']].corr()
		print '\n'	
		print self.data[['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']].corr()
		print '\n'

	def cell_type_corr(self):
		'''Compute correlation coefficient between cell type'''
		print "The correlation matrix for cell types is:\n"
		celltype = DataFrame({'HL60':self.data['HL60_0_hrs'].append(self.data['HL60_24_hrs']),'U937':self.data['U937_0_hrs'].append(self.data['U937_24_hrs']),'Jurkat':self.data['Jurkat_0_hrs'].append(self.data['Jurkat_24_hrs call']),'NB4':self.data['NB4_0_hrs'].append(self.data['NB4_24_hrs'])})
		print celltype.corr()
		print '\n'

	def select_gene(self):
		'''selcet gene change little in smaples'''
		print "The selected genes are:\n"
		self.data['std'] = self.data[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs','U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs','Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs call','NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']].std(axis=1)
		data_sort = self.data.sort('std')
		print data_sort.head(n=10)['Gene Accession Number']
		print '\n'

	def QE(self):
		print "The genes are:\n"
		data_select = self.data[(self.data['NB4_24_hrs'] > 2*self.data['NB4_0_hrs']) & (self.data['HL60_24_hrs'] > 2*self.data['HL60_0_hrs']) & (self.data['U937_24_hrs'] > 2*self.data['U937_0_hrs']) & (self.data['Jurkat_24_hrs call'] > 2*self.data['Jurkat_0_hrs'])]
		for index, row in data_select.iterrows():
			print row['Gene Accession Number']
		print '\n'	

	def QF(self):
		print "The differentially regulated genes are:\n"
		data_select = self.data[(self.data['HL60_0_hrs'] > 2*self.data['U937_0_hrs']) | (2*self.data['HL60_0_hrs'] < self.data['U937_0_hrs'])]                
		for index, row in data_select.iterrows():
                        print row['Gene Accession Number']
                print '\n'

if __name__ == '__main__':
	ea = ExpressionAnalyzer()
	ea.import_file('../data_set_HL60_U937_NB4_Jurkat.txt')
	ea.count_gene()
	ea.time_corr()
	ea.cell_type_corr()
	ea.select_gene()
	ea.QE()
	ea.QF()
