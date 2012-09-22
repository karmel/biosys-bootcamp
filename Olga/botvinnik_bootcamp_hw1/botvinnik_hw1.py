#/usr/bin/env python

# R-like DataFrames
import pandas

# Parse command line arguments
import argparse

# Regular expressions
import re

# Matlab-like matrix operations
import numpy as np

'''
UCSD Bioinformatics and Systems Biology Bootcamp
Homework 1

Author: Olga Botvinnik

The purpose of this program is to import a dataset using the pandas package, 
convert it into a DataFrame, and query properties of the data, such as:
1. How many distinct genes are represented?
2. Which two time points are the most highly correlated for each 
   sample type?
3. Which two cell types are the most similar?
4. What are top 10 genes with the least amount of variation? (used for 
   normalization)
5. Do any genes show two-fold higher expression at time point X vs time 
   point Y, for all cell types? If so, which?
6. Which genes are differentially regulated in sample A as compared to 
   sample B, at time point X?

Example run:
python botvinnik_hw1.py -f data/data_set_HL60_U937_NB4_Jurkat.csv -N 10 -t0 0_hrs -t1 24_hrs --sampleA HL60 --sampleB U937 --timeX _0_hrs
'''

class CommandLine:
	def __init__(self, inOpts=None):
		self.parser = argparse.ArgumentParser(description = '''Import a dataset
			and query properties such as number of distinct genes, the two
			most highly correlated time points, which cell types are most
			similar, top N genes with least amount of variation, any genes
			with 2-fold higher expression at time point X vs time point Y,
			for all cell types, differentially expressed genes in sample A vs
			sample B at time point X.''',
			add_help=True, prefix_chars='-', 
			usage='%(prog)s --file FILE [options]')
		self.parser.add_argument('--file', '-f', action='store',
			help='the datafile to import and query', required=True)
		self.parser.add_argument('--top-N-least-var', '-N', action='store',
			help='state the top N genes with the least amount of variation.')
		self.parser.add_argument('--time0', '-t0', action='store',
			help='For finding two-fold higher expression for all samples,\
			the first time point')
		self.parser.add_argument('--time1', '-t1', action='store',
			help='For finding two-fold higher expression for all samples,\
			 the second time point')
		self.parser.add_argument('--sampleA', '-sA', action='store',
			help='For finding differentially expressed genes between sample A\
			and sample B at time point X')
		self.parser.add_argument('--sampleB', '-sB', action='store',
			help='For finding differentially expressed genes between sample A\
			and sample B at time point X')
		self.parser.add_argument('--timeX', '-tX', action='store',
			help='For finding differentially expressed genes between sample A\
			and sample B at time point X')
		if inOpts is None:
			self.args = vars(self.parser.parse_args())
		else:
			self.args = vars(self.parser.parse_args(inOpts))

	def do_usage_and_die (self, str):
		'''
		If a critical error is encountered, where it is suspected that the 
		program is not being called with consistent parameters or data, this
		method will write out an error string (str), then terminate execution 
		of the program.
		'''
		import sys
		print >>sys.stderr, str
		self.parser.print_usage()
		return 2

class Usage(Exception):
	'''
	Used to signal a Usage error, evoking a usage statement and eventual
	exit when raised
	'''
	def __init__(self, msg):
		self.msg = msg

class ExpressionAnalyzer(object):
	def __init__(self, cl_args):
		'''
		Given a dictionary of command line arguments ("cl_args"), initialize the 
		ExpressionAnalyzer object
		'''
		# --- Creation and modification of microarray counts DataFrame --- #
		# Import the file and turn it into a pandas DataFrame object
		self.import_file(cl_args['file'])
		# print self.df

		self.genbank = self.df['Gene Accession Number']

		# Remove columns that have 'object' as the data type 
		# (should be a number type, in this case, all the 'call' columns
		# and the 'Gene Accession Number' columns are ignored)
		self.df = self.df[self.df.columns[~(self.df.dtypes == 'object')]]
		
		# The original values for this particular dataset were integers,
		# ans we need to conver them into floats to calculate the correlation
		# print df_clean
		self.df = self.df.applymap(float)

		# --- Creation of correlation DataFrame --- #
		# Get a matrix of correlations
		self.correl = self.df.corr()

		# Find the samples that are maximally correlated within each sample type
		sample_types = self.find_max_correl_samples()

		# ---  Saving other command line variables --- #

		#######################################################################
		# 1. QUESTION: is there a cute single-line way to do this? In R I would do
		# something along the lines of:
		# N = ifelse(cl_args['top_N_least_var']==None, None, cl_args['top_N_least_var'])
		#######################################################################		
		if cl_args['top_N_least_var'] is not None:
			self.N = int(cl_args['top_N_least_var'])
		else:
			self.N = None

		# Compare all samples of one type at two different time points. Which genes
		# are differentially expressed for ALL sample types at these time points?
		self.t0 = cl_args['time0']
		self.t1 = cl_args['time1']

		# Compare sample A and sample B at time X to find 
		# differentially expressed genes
		self.sampleA = cl_args['sampleA']
		self.sampleB = cl_args['sampleB']	
		self.timeX = cl_args['timeX']

	def import_file(self, filename):
		'''
		Imports a data file and returns a pandas DataFrame
		'''
		self.df = pandas.read_csv(filename, index_col=0)

	def print_num_distinct_genes(self):
		'''
		Given a DataFrame, use the DataFrame.index.unique() function to print
		the number of distinct genes
		'''
		print 'There are', len(self.df.index.unique()), 'distinct genes'

	def df_remove_re_match_col(self, df, s):
		'''
		Given a dataframe, remove the columns that match this string pattern,
		which will be converted into a regular expression object.
		Return the modified DataFrame.
		'''
		# create the regular expression function to get rid of the pattern s
		r = re.compile(s)
		vmatch = np.vectorize(lambda x:bool(r.match(x)))
		return df[df.columns[~vmatch(df.columns)]]

	def df_only_re_search_col(self, df, s):
		'''
		Given a DataFrame, take only the columns that can be re.search'd with
		this string pattern, which will be converted into a regular expression object.
		Return the modified DataFrame.
		'''
		r = re.compile(s)
		vsearch = np.vectorize(lambda x:bool(r.search(x)))
		return df[df.columns[vsearch(df.columns)]]

	def df_only_re_search_row(self, df, s):
		'''
		Given a DataFrame, take only the columns that can be re.search'd with
		this string pattern, which will be converted into a regular expression object.
		Return the modified DataFrame.
		'''
		r = re.compile(s)
		vsearch = np.vectorize(lambda x:bool(r.search(x)))
		return df.ix[vsearch(df.index)]

	def find_max_correl_samples(self):
		'''
		Given a data frame, find the two samples that have the highest pearson
		correlation. Also, set unique sample types variable 'self.sample_types', 
		without the time point trailing information, e.g. without '_0_hrs'
		'''
		# Now need to get correlations within cell types only.
		#### Note: should add command line flag to specify '_0_hrs' or other ####
		#### common sample name endings to get the core sample names ####
		r = re.compile('(.+)_0_hrs')
		vsearch = np.vectorize(lambda x:bool(r.search(x)))

		# format() forces the pandas Index format into a list of strings
		# which can then be easily searched
		sample_types = self.correl.columns[vsearch(self.correl.columns.format())]
		# print sample_types

		# Create a new array to store the 'cleaned' sample type names, without
		# the trailing '_0_hrs'
		sample_types_clean = []
		for s in sample_types:
			m = re.match('(.+)_0_hrs', s)
			# if m is not None:
			sample_types_clean.append(m.group(1))
		sample_types = sample_types_clean

		for s in sample_types:
			# Take a subset of the correlation data frame that only includes
			# the correlations within these samples
			# Note: for some reason, DataFrame.lookup wasn't working for me
			#       to obtain specific row,column pairs.

			# Get just the columns with this sample type
			correl_mini = self.df_only_re_search_col(self.correl, s)
			# print correl_mini

			# Get just the rows with this sample type
			correl_mini = self.df_only_re_search_row(self.correl, s)
			# print correl_mini

			# We want the samples with the maximum correlation, but any sample
			# will be maximally correlated with itself, so lets set all the
			# correl = 1 to 0:
			correl_mini = correl_mini.replace(1, 0)

			max_value = 0
			max_index = None
			max_column = None
			# print correl_mini
			for c in correl_mini.columns:
				# print 'c', c
				mx = correl_mini[c].max()
				# print 'mx', mx
				if mx > max_value:
					max_value = mx
					max_index = correl_mini[c].argmax()
					max_column = c
					max_row = correl_mini.index[max_index]
			print 'The most correlated samples within %s are: %s and %s (correl = %.4f).'\
			% (s, max_column, max_row, max_value)
		self.sample_types = sample_types

	def find_similar_cell_types(self):
		'''
		Given a DataFrame, find the two samples that consistently have the highest
		pearson correlation.
		'''
		max_sumsq_correl = 0
		max_pair_correl = (None, None)

		# Create a list of unordered sets to check against to make sure
		# we haven't already made this comparison
		pairs_already_explored = []

		for s1 in self.sample_types:
			for s2 in self.sample_types:
				# If we've already explored this pair or the pair consists
				# of the same sample types
				pair = set([s1, s2])
				if s1 == s2 or pair in pairs_already_explored:
					continue
				# print pair
				correl_mini = self.df_only_re_search_col(self.correl, s1)

				# Feed this new correl_mini dataframe into the function
				# to filter rows
				correl_mini = self.df_only_re_search_row(correl_mini, s2)
				# print correl_mini

				# Square each element of the array to more strictly penalize
				# float distance from 1
				sumsq_correl = np.square(correl_mini).sum().sum()			
				if sumsq_correl > max_sumsq_correl:
					max_sumsq_correl =  sumsq_correl
					max_pair_correl = pair
				
				# Mark this pair as already explored
				pairs_already_explored.append(pair)

		max_pair_correl = list(max_pair_correl)
		print 'The sample types that are most highly correlated are %s and %s (sum-squared of correl matrix: %.4f).'\
		% (max_pair_correl[0], max_pair_correl[1], max_sumsq_correl)

	def find_least_variant_genes(self):
		'''
		Given a DataFrame and an integer N, find the top N genes with the least
		variance
		'''
		# Calculate variances on a per-gene basis, across all samples
		variances = self.df.var(1)

		# Get the names of the top N genes with the least variance without
		# touching the original order
		least_variant = variances.order()[range(self.N)]

		print 'The least variant genes are:'
		for i in range(self.N):
			print '\t%s\t(variance = %.4f)' % (least_variant.index[i], least_variant[i])

	def find_two_fold_diff_expr_all(self):
		'''
		Given a DataFrame and two time points t0 and t1, find the union of genes thare
		are differentally expressed in all sample types
		'''
		genes_2fold_diff_expr = None
		for s in self.sample_types:
			# Get only the columns corresponding to this sample type and corresponding
			# time points
			df_mini = self.df_only_re_search_col(self.df, s)
			df_t0 = self.df_only_re_search_col(df_mini, self.t0)
			df_t1 = self.df_only_re_search_col(df_mini, self.t1)

			# Get the column names for indexing
			t0_name = df_t0.columns[0]
			t1_name = df_t1.columns[0]

			# print t0_name, t1_name

			#######################################################################
			# 2. QUESTION: Concatenating these dataframes created NaN's. How can you
			# "column-bind" two dataframes in pandas? Concatenate tries too hard to
			# figure out how to match up the indices, and I want to just force
			# column binding without paying attention to the row names, since in
			# this case I know they're exactly the same.
			#######################################################################			

			# Concatenating these dataframes created NaN's so I'm sticking
			# to just analyzing df_mini
			s_genes_2fold = df_mini[df_mini[t1_name] >= 2*df_mini[t0_name]].index
			
			if genes_2fold_diff_expr is None:
				# Create a set, which is an unordered collection of unique objects
				genes_2fold_diff_expr = set(s_genes_2fold)
			else:
				# Intersect our old set with the new set to get rid of genes
				# that are not common to both
				genes_2fold_diff_expr = genes_2fold_diff_expr.intersection(set(s_genes_2fold))

			# print 'len(s_genes_2fold)', len(s_genes_2fold)
			# print 'len(genes_2fold_diff_expr)', len(genes_2fold_diff_expr)
			# print 'len(genes_2fold_diff_expr & s_genes_2fold)', len(genes_2fold_diff_expr & s_genes_2fold)
		print 'There are %d genes that are two-fold differentially expressed in all sample types.'\
		% (len(genes_2fold_diff_expr))

	def diff_expr_sampleA_sampleB_timeX(self):
		'''
		Given the DataFrame of expression data, check which genes are 2-fold
		or more differentially expressed in sample A over sample B.
		'''

		# Make a 'mini' dataframe with just the samples of this time point
		df_mini = self.df_only_re_search_col(self.df, self.timeX)

		# Make dataframes with just
		df_A = self.df_only_re_search_col(df_mini, self.sampleA)
		df_B = self.df_only_re_search_col(df_mini, self.sampleB)

		A_name = df_A.columns[0]
		B_name = df_B.columns[0]

		diff_expr_genes = self.genbank[df_mini[A_name] >= 2*df_mini[B_name]]
		print diff_expr_genes.head()

		txt_file = 'diff_expr_%s_vs_%s_at_t=%s.txt' % \
			(self.sampleA, self.sampleB, self.timeX)
		diff_expr_genes.tofile(txt_file, sep="\n")

		print 'There are %d differentially expressed genes between %s vs %s at time %s.'\
		% (len(diff_expr_genes), self.sampleA, self.sampleB, self.timeX)
		print 'GenBank accession codes written to: %s' % txt_file

def import_and_query():
	cl = CommandLine()
	try:
		# Debugging statement: print out the command line arguments
		# to make sure they were parsed correctly
		print cl.args

		# Initialize ExpressionAnalyzer object (includes converting the data file 
		# into a data frame)
		ea = ExpressionAnalyzer(cl.args)

		# --- A. How many distinct genes? --- #
		ea.print_num_distinct_genes()

		# --- B. Which two time points are most highly correlated --- #
		# --- for each cell type?                                 --- #
		ea.find_max_correl_samples()

		# --- C. Which two cell types are the most similar? --- #
		ea.find_similar_cell_types()

		# --- D. Which are the top N genes with the lowest variation? --- #
		if ea.N is not None:
			ea.find_least_variant_genes()

		# --- E. Do any genes show two-fold higher expression at time1 vs  --- #
		# --- vs time0 hours for ALL cell types? (here, time1=24, time0=0) --- #
		if ea.t0 is not None and ea.t1 is not None:
			ea.find_two_fold_diff_expr_all()

		# --- F. Which genes are differentially expressed in sampleA (HL60) --- #
		# --- vs sampleB (U937) at timeX (0 hours)?                         --- #
		if ea.sampleA is not None and ea.sampleB is not None and ea.timeX is not None:
			ea.diff_expr_sampleA_sampleB_timeX()

	# If not all the correct arguments are given, break the program and
	# show the usage information
	except Usage, err:
		cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
	import_and_query()