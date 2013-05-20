#################################################################
########################## ANSWERS ##############################
#################################################################
# Bootcamp Day 1 Homework | Nathan Mih
# A: There are 5001 distinct genes represented in the data set.
# B: For each cell type, the highest correlations are between the following time points: 
#		HL60: 0 and 4, U937: 0.5 and 4, NB4: 24 and 48, Jurkat: 0.5 and 4
# C: The two cell types that are most similar at the 0 time point are Jurkat and U937.
# D: The 10 genes with the least variance: 
# Gene Accession Number
# H85478                   2.735294
# H88261                   2.882353
# M95586_r_i               3.632353
# X73424_i                 3.779412
# R94219                   4.183824
# U39231                   4.235294
# H82137                   4.264706
# L20814                   4.279412
# R39006                   4.360294
# X86401                   4.404412
# E: The following genes show two-fold higher expression: 
# 		H56092, H63361, H68366, H88978, L14075, M18700, M60828, R35903, R49046, R49688, R55782_i, R67283, T51288, U12978, U18934_i, V00563, X00090_i, X06700, X51420, X58987, Z22534_i, Z48051
# F: There are over 4000 genes that are differentially regulated. See file different.csv for all the genes.
###################################################################
# G: According to the DAVID ontology analyzer, there are 1775 DAVID IDs from the list. 
###################################################################

# Well-commented. Props for that!

from __future__ import division # Isn't there something rewarding about importing from THE FUTURE?
from pandas import *
import numpy as np
import pandas as pd


class ExpressionAnalyzer(object):
	'''
	Responsible for importing and analyzing gene expression data.
	Expects data with named samples as columns and genes as rows.
	'''

	def import_file(self, file):
		''' Import data from file and return a DataFrame object. '''
		# read_csv will read entire file and parse it to a DataFrame
		# header = 0 because first row contains headers for columns
		# index_col = 1 for setting 'Gene Accession Number' to define rows
		self.df = pd.read_csv(file, header = 0, index_col = 1)
		return self.df

	def num_genes(self, column):
		''' Tells us how many distinct genes are represented in the data set. '''
		# get rid of duplicates in 'Gene Description', return length
		no_dupe = set(self.df[column])
		return len(no_dupe)

	def time_corr(self):
		''' Which two time points are the most highly correlated for each cell type? '''
		'''	Which two cell types are the most similar? '''
		# creates a file with an array of correlations
		# need to implement way to read and compare correlations (not currently here)
		correlations = self.df.corr(method = 'pearson') # Parenthesis around self.df not necessary here or elsewhere
		# outputs correlations to a csv file
		correlations.to_csv('hw1_correlations.csv')

	def least_variance(self):
		'''	Which 10 genes change very little across samples? '''
		# look at variances for all genes
		variances = (self.df).var(axis = 1)
		# sort and return the 10 with the least variance
		variances.sort()
		return variances[:10]

	# not very general code here
	# At least you recognize it :) Would you know how to generalize it the next time around?
	def two_fold(self, hl0, hl24, u0, u24, nb0, nb24, j0, j24):
		'''	Do any genes show two-fold higher expression at 24 hours versus 0 hours for all four cell types?  '''
	 	'''	Which genes are differentially regulated (at least two-fold higher or lower) in HL60 cells as compared to U937 cells at 0 hours? '''
		# again, need to implement way to read array instead of manually looking
		# idea: loop thru all genes - if # at 24 is greater than 2x # at 0 for all cell types then input gene name into a list
		# twofold will create an array with boolean values
		twofold = (abs((self.df)[hl0]/(self.df)[hl24]) > 2) & (abs((self.df)[u0]/(self.df)[u24]) > 2) & (abs((self.df)[nb0]/(self.df)[nb24]) > 2) & (abs((self.df)[j0]/(self.df)[j24]) > 2)
		twofold.sort()
		# outputs array with true values meaning they show higher expression
		twofold.to_csv('hw1_twofold.csv')

		# same principle as above array
		different = (abs((self.df)[hl0]/(self.df)[u0]) > 2) | (abs((self.df)[u0]/(self.df)[hl0]) > 2)
		different.sort()
		# outputs array with true values meaning that they are differentially regulated
		different.to_csv('hw1_different.csv')


# Add if __name__ == '__main__': and indent. Google for why that is important.
# answer output section
ea = ExpressionAnalyzer()
ea.import_file('data_set_HL60_U937_NB4_Jurkat.csv')

print "Bootcamp Day 1 Homework | Nathan Mih"

# Question A
print "A: There are %r distinct genes represented in the data set." % (ea.num_genes('Gene Description'))

# Question B
ea.time_corr()
# looked at correlations manually and recorded the highest ones between cell types
# HL60: 0 and 4
# U937: 0.5 and 4
# NB4: 24 and NB4_48_hrs
# Jurkat: 0.5 and 4
print "B: For each cell type, the highest correlations are between the following time points: HL60: 0 and 4, U937: 0.5 and 4, NB4: 24 and 48, Jurkat: 0.5 and 4"

# Question C
# looked at correlations again and looked for highest correlations at 0 hours
# Jurkat and U937 are most similar (highly correlated)
print "C: The two cell types that are most similar at the 0 time point are Jurkat and U937."

# Question D
print "D: The 10 genes with the least variance: \n%r" % ea.least_variance()

# Question E
ea.two_fold('HL60_24_hrs', 'HL60_0_hrs', 'U937_0_hrs', 'U937_24_hrs', 'NB4_0_hrs', 'NB4_24_hrs', 'Jurkat_0_hrs', 'Jurkat_24_hrs')
print "E: The following genes show two-fold higher expression: H56092, H63361, H68366, H88978, L14075, M18700, M60828, R35903, R49046, R49688, R55782_i, R67283, T51288, U12978, U18934_i, V00563, X00090_i, X06700, X51420, X58987, Z22534_i, Z48051"

# Question F
print "F: There are over 4000 genes that are differentially regulated. See file different.csv for all the genes."