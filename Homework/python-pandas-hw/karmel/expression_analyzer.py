from __future__ import division
import sys
from pandas import DataFrame
from pandas.io import parsers
import numpy as np

class ExpressionAnalyzer(object):
	'''
	Responsible for importing and analyzing gene expression data.
	Expects data with named samples as columns and genes as rows.
	'''
	
	def import_file(self, filename, sep='\t', header=0):
		''' Import data from file and return a DataFrame object. '''
		return parsers.read_csv(filename, sep=sep, header=header)
		
	def count_distinct_genes(self, data):
		''' Returns the number of distinct gene names in the dataset. '''
		return len(set(data['Gene Description']))
		
	def corr_timepoints(self, data, cell_type='', timepoints=None):
		''' Determines the correlation of the passed timepoint and cell type. '''
		key_set = ['{}_{}_hrs'.format(cell_type, tp) for tp in timepoints]
		return np.array(data[key_set].corr())

	def corr_cells(self, data, cell_types=None):
		''' Determines the correlation of the passed cell types and timepoints. 
		cell_types should be a tuple of cell, [timepoints] tuples.
		Strategy: subselect the relevant columns for each cell type, flatten into a 1D array,
		and get the correlation coefficient of each pair of 1D cell_type arrays.
		'''
		# First, prep 1D vectors for each cell type
		subsets = []
		for cell_type, timepoints in cell_types:
			key_set = ['{}_{}_hrs'.format(cell_type, tp) for tp in timepoints]
			subsets.append(data[key_set].unstack())
		
		# Pairwise correlation of the flattened, 1D cell vectors
		correls = np.corrcoef(subsets)
		row, col = self.get_highest_correlation(correls)
		
		# Given the highest correlation, return the cell type names
		return cell_types[row][0], cell_types[col][0]
		
	def get_highest_correlation(self, correls):
		''' Given a 2D array of correlation coefficients, 
			return the x,y index of the highest value.
		'''
		# np.triu gives the upper-triangle of an array; 
		# k=1 shifts the selection so that we don't get the diagonal
		highest = np.triu(correls, k=1).argmax()
		row = int(highest/correls.shape[1])
		col = highest % correls.shape[0]
		return row, col

	def get_var(self, data, cell_types):
		''' Returns 1D array of variation per gene.'''
		return data[self.get_key_set(cell_types)].var(axis=1)
	
	def genes_with_min_change(self, data, cell_types, k=10):
		''' Return gene list sorted by variation.'''
		data['var'] = self.get_var(data, cell_types)
		return data.sort_index(by='var')[['Gene Accession Number','Gene Description']][:k]
	
	def count_genes_higher(self, data, cell_type_names, t1, t2, fold_change=2):
		''' Return count of genes fold_change higher in t1 than t2 for all cell_type_names.'''
		higher = None
		for cell in cell_type_names:
			higher_bool = data['{}_{}_hrs'.format(cell, t1)] > fold_change*data['{}_{}_hrs'.format(cell, t2)]
			higher_set = set(data[higher_bool]['Gene Accession Number'])
			if higher is None: higher = higher_set
			else: higher = higher_set & higher
		return len(higher)
	
	def genes_differential(self, data, cell1, cell2, timepoint, fold_change=2):
		''' Return list of genes fold_change higher between cell1 and cell2 at given timepoint. '''
		higher_bool = data['{}_{}_hrs'.format(cell1, timepoint)] > fold_change*data['{}_{}_hrs'.format(cell2, timepoint)]
		lower_bool = data['{}_{}_hrs'.format(cell1, timepoint)]*fold_change < data['{}_{}_hrs'.format(cell2, timepoint)]
		# Return accession numbers less the suffix after the underscore, which is some sort of isoform id
		return map(lambda x: x.split('_')[0], data[higher_bool & lower_bool]['Gene Accession Number'])
		
	key_set = None
	def get_key_set(self, cell_types):
		''' Get full set of numeric column names.'''
		if not self.key_set:
			self.key_set = ['{}_{}_hrs'.format(cell_type, tp) 
								for cell_type, timepoints in cell_types
								for tp in timepoints]
		return self.key_set
		
class BroadExpressionAnalyzer(ExpressionAnalyzer):
	cell_types = (('HL60', [0, .5, 4, 24]),
				('U937', [0, .5, 4, 24]),
				('NB4', [0, 5.5, 24, 48]),
				('Jurkat', [0, .5, 4, 24],))
	
if __name__ == '__main__':

	try: filename = sys.argv[1]
	except IndexError: filename = 'data_set_HL60_U937_NB4_Jurkat.txt'
	ea = BroadExpressionAnalyzer()

	data = ea.import_file(filename)
	print 'Part A: Distinct Genes'
	print 'There are {} distinct genes.'.format(ea.count_distinct_genes(data))

	print 'Part B: Correlated Timepoints'
	for cell, timepoints in ea.cell_types:
		correls = ea.corr_timepoints(data, cell, timepoints)
		row, col = ea.get_highest_correlation(correls)
		print 'For cell type {}, timepoints {} and {} are most correlated ({})'.format(
									cell, timepoints[row], timepoints[col], 
									correls[row][col])
	
	print 'Part C: '
	print '{} and {} are the two most correlated cell types.'.format(*ea.corr_cells(data, ea.cell_types))
	
	print 'Part D: Stable Genes'
	print 'Genes with minimal change include:'
	for _, row in ea.genes_with_min_change(data, ea.cell_types).iterrows(): 
		print '\t{}:{}'.format(row['Gene Accession Number'],row['Gene Description'],)
	
	print 'Part E: Count of Two-fold Higher Genes'
	print '{} genes are 2x greater in 24 hours than 0 hours in all cell types.'.format(
				ea.count_genes_higher(data, zip(*ea.cell_types)[0], 24, 0))
				
	print 'Part F: Differing Gene List'
	print 'Genes that are 2x different in HL60 and U937 cells at 0 hours:'
	print ', '.join(ea.genes_differential(data, 'HL60', 'U937', 0))
	
	
