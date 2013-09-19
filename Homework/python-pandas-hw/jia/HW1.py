

from pandas import *
from numpy import *

class ExpressionAnalyzer(object):
	'''
	Responsible for importing and analyzing gene expression data.
	Expects data with named samples as columns and genes as rows.
	'''
    def __init__(self):
	'''initiation'''
        self.fileadd = {}

    def import_file(self, filename):
	'''import file'''
        self.fileadd = pandas.read_csv(filename, sep='\t', index_col=False)
    
    def count_gene_number(self):
	'''count gene number'''
        s=Series(self.fileadd['Gene Accession Number'])
        print(len(s.unique()))
    
    def cal_cor(self):
	'''calculate the correlation'''
        print (self.fileadd[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs']].corr(),'\n')
        print (self.fileadd[['U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs']].corr(),'\n')
	print(self.fileadd[['NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs']].corr(),'\n')
        print(self.fileadd[['Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrscall']].corr(),'\n')

    def cellline_cor(self):
	'''calculate cell lines correlation'''
	df=DataFrame({'HL60':self.fileadd['HL60_0_hrs'].append(self.fileadd['HL60_24_hrs']),'U937':self.fileadd['U937_0_hrs'].append(self.fileadd['U937_24_hrs']),'NB4':self.fileadd['NB4_0_hrs'].append(self.fileadd['NB4_24_hrs']),'Jurkat':self.fileadd['Jurkat_0_hrs'].append(self.fileadd['Jurkat_24_hrscall'])})
	print(df.corr(),'\n')

	
    def ten(self):
	self.fileadd['std'] = self.fileadd[['HL60_0_hrs','HL60_0.5_hrs','HL60_4_hrs','HL60_24_hrs','U937_0_hrs','U937_0.5_hrs','U937_4_hrs','U937_24_hrs','Jurkat_0_hrs','Jurkat_0.5_hrs','Jurkat_4_hrs','Jurkat_24_hrs call','NB4_0_hrs','NB4_5.5_hrs','NB4_24_hrs','NB4_48_hrs','NB4_72_hrs']].std(axis=1)
	sort = self.fileadd.sort('std')
	print(sort.head(n=10)['Gene Accession Number'],'\n')
    

    def find_two_folders(self):
        for i in range(len(fileadd['Gene Description'])):
            if self.fileadd['HL60_0_hrs'][i]>2*self.fileadd[HL60_24_hrs] and self.fileadd[U937_0_hrs][i]>2*self.fileadd[U937_24_hrs] and self.fileadd[NB4_0_hrs][i]>2*self.fileadd[NB4_24_hrs] and self.fileadd[Jurkat_0_hrs][i]>2*self.fileadd[Jurkat_24_hrscall]:
                print(self.fileadd['Gene Accession Number'][i],'\n')
    

    def diff_regular(self):
	'''find different regulated gene'''
        for i in range(len(self.fileadd['Gene Description'])):
            if self.fileadd['HL60_0_hrs'][i]>2.0*self.fileadd['U937_0_hrs'][i] or self.fileadd['HL60_0_hrs'][i]<0.5*self.fileadd['U937_0_hrs']:
                print(self.fileadd['Gene Accession Number'][i],'\n')
                      

ea=ExpressionAnalyzer()
ea.import_file('./data4.csv')
ea.count_gene_number()
ea.cal_cor()
ea.cellline_cor()
ea.ten()
ea.find_two_folders()
ea.diff_regular()





