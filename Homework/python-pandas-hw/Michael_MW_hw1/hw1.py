from __future__ import division
import pandas as pd
import numpy as np
from pandas import *

class ExpressionAnalyzer:
# read file
	def read_csv(self, filename):
		self.dataset = pd.read_csv(filename, sep='\t')
		
#A. How many distinct genes are represented in the data set?
	def Q_A(self):
		self.gene_set = len(set(self.dataset['Gene Description']))
        	print self.gene_set		
        	
#B. Which two time points are the most highly correlated for each cell type? 
	def Q_B(self):
		print self.dataset[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs']].corr()
		print self.dataset[['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs']].corr()
		print self.dataset[['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']].corr()
		print self.dataset[['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs']].corr()
		

#C. Which two cell types are the most similar?
	def Q_C(self):	
		df=DataFrame({'HL60':self.dataset['HL60_0_hrs'].append(self.dataset['HL60_24_hrs']),
		              'U937':self.dataset['U937_0_hrs'].append(self.dataset['U937_24_hrs']),
		              'NB4':self.dataset['NB4_0_hrs'].append(self.dataset['NB4_24_hrs']),
		              'Jurkat':self.dataset['Jurkat_0_hrs'].append(self.dataset['Jurkat_24_hrs'])})
		print df.corr()					

#D.It is often useful to know which genes change very little across samples for the sake of normalization or calibration. Based on this data set, what are ten good candidates for genes to use to calibrate machinery or analyses across all these samples?
	def Q_D(self):
		self.dataset['std'] = self.dataset[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs',
		'U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs',
		'Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs',
		'NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']].std(axis=1)
		sort = self.dataset.sort('std')
		candidates=10
		print sort.head(candidates)['Gene Accession Number']
		
#E. Do any genes show two-fold higher expression at 24 hours versus 0 hours for all four cell types? If so, which ones?
	def Q_E(self):
		self.cl_list = ['HL60','U937','NB4','Jurkat']
		self.count =0
		for gene in range(len(self.dataset['Gene Description'])):
			self.mr_true = 0
			for cl in self.cl_list:	
				if self.dataset[cl+'_24_hrs'][gene] > 2.0*self.dataset[cl+'_0_hrs'][gene]:
					self.mr_true = 1
				else: 
					self.mr_true = 0
				self.mr_true *= self.mr_true
			if self.mr_true == 1:
				self.count += 1
		print self.count
	
#F. Which genes are differentially regulated (at least two-fold higher or lower) in HL60 cells as compared to U937 cells at 0 hours? 
	def Q_F(self):
		filehandle= open ('genelist.txt','w')
		fold_change =1.0
        	for gene in range(len(self.dataset['Gene Accession Number'])):
        		fold_change = self.dataset['HL60_0_hrs'][gene] / self.dataset['U937_0_hrs'][gene]
        		if (fold_change > 2.0 and self.dataset['HL60_0_hrs'][gene]>0.0) or (fold_change < 0.5 and fold_change >= 0.0 and self.dataset['HL60_0_hrs'][gene]>0.0):
                		filehandle.write(self.dataset['Gene Accession Number'][gene]+'\n')
	


#G.Take the list of Gene Accession codes from (F), and run them through the DAVID ontology analyzer. (at http://david.abcc.ncifcrf.gov/summary.jsp . These are GenBank Accession codes.) Are there any enriched ontology terms?
	#see ans.txt

#main
def main():
	
	#configuration:
	filename='data_set_HL60_U937_NB4_Jurkat.csv'
	hw1 = ExpressionAnalyzer()
	hw1.read_csv(filename)

	#call solution:
	'''	hw1.Q_A()
	hw1.Q_B()
	hw1.Q_C()
	hw1.Q_D()
	hw1.Q_E()'''
	hw1.Q_F()
	
	''' #for some reason the following is not working:
        for problem in 'ABCDEF':
        	print 'Q_' + problem
        	hw1.__getattribute__('Q_'+ problem)()
	'''
	
if __name__ == '__main__':
	main()

	
