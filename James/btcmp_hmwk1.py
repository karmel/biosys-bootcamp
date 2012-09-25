import pandas as pd
import numpy as np
import sys
import re
import math
# Object oriented goodness. Awesome!
class ExpressionAnalyzer(object):
	"""
	Object for analyzing gene expression files of the same general format as
	http://www.broadinstitute.org/mpr/publications/projects/SOM_Methods_and_Applications/data_set_HL60_U937_NB4_Jurkat.txt
	"""
	
	def import_file(self, exp_file):
		"""
		Imports csv file (could be modified to detect and handle other file types).
		File contents are stored in 'data' attribute.
		Also reads first line into a separate 'header' attribute. #This is probably unnecessary
		"""
		extension = re.search('\.([a-z]+)$', expression_file)
		extension = 'csv'
		self.delimiter = ''
		if extension is 'txt':
			self.delimiter = '\t'
		elif extension is 'csv':
			self.delimiter = ','
		self.data = pd.read_csv(exp_file, sep=self.delimiter, header=0)
		with open(exp_file, 'r') as f:
			 self.header = f.readline()
		 
	def count_genes(self):
		"""
		Parameters: none
		Returns: the number of genes (assuming each record in the file is a unique gene)
		"""
		# Not all genes were unique in the file --> must filter out dupes.
		return len(self.data)
		# when might this be insufficient?
		
	def field_index(self, field_name):
		"""
		Parameters: name of a column
		Returns: index of column with a given name. If not found, returns -1
		"""
		# Great commenting, with params and return notes on every object. Love it.
		result_index = -1
		pattern = "(.*?)"
		pattern += self.delimiter
		all_column_names = re.findall(pattern, self.header)
		# Check out the pandas built in .filter() method too for regex
		for index in range(0,len(all_column_names)):
			if all_column_names[index] == field_name:
				result_index = index
				break
		return result_index
		
	def relevant_column_names(self, prefix):
		"""
		Parameters: a cell type prefix
		Returns: a list of columns containing the specified prefix
		"""
		# Awesome. This is a method I should have factored out in my solution,
		# but didn't. Nice that you did.
		pattern = ",("
		pattern += prefix
		pattern += ".*?)["
		pattern += self.delimiter
		pattern += "|\n]"
		names = re.findall(pattern, self.header)
		return names
		
	def most_correlated_time_points(self, list_of_prefixes):
		"""
		Parameters: a list of cell type prefixes.
		Returns: a list of tuples. Each tuple holds a cell type prefix
			and the two time points for which that cell type's expression values are most correlated.
		"""
		time_point_pairs = []
		for prefix in list_of_prefixes:
			max_correlation = 0
			time_point_pair = ["", ""]
			col_names = self.relevant_column_names(prefix)
			number_of_columns = len(col_names)
			for num in range(0,number_of_columns):
				for other_num in range(num + 1, number_of_columns):
					column_name_1 = col_names[num]
					column_name_2 = col_names[other_num]
					col_1 = self.data[column_name_1]
					col_2 = self.data[column_name_2]
					correlation = col_1.corr(col_2)
					if correlation > max_correlation:
						max_correlation = correlation
						time_1 = re.findall('_([0-9]+\.?[0-9]?)_', column_name_1)
						time_2 = re.findall('_([0-9]+\.?[0-9]?)_', column_name_2)
						time_point_pair = [time_1, time_2]
			time_point_pairs.append([prefix, time_point_pair])
		return time_point_pairs

	def most_similar_cell_types(self, list_of_prefixes):
		"""
		Parameters: a list of cell type prefixes
		Returns: two prefixes corresponding to the two cell types with 
			highest correlation between expression values, averaged across all common time points
		"""
		highest_mean_corr = 0
		cell_types = ["", ""]
		number_of_prefixes = len(list_of_prefixes)
		for p_num in range(0, number_of_prefixes):
			for other_p_num in range(p_num + 1, number_of_prefixes):
				pre_1 = list_of_prefixes[p_num]
				pre_2 = list_of_prefixes[other_p_num]
				col_names_1 = self.relevant_column_names(pre_1)
				col_names_2 = self.relevant_column_names(pre_2)
				correlations = []
				for col_name_1 in col_names_1:
					for col_name_2 in col_names_2:
						time_1 = re.findall('_([0-9]+\.?[0-9]?)_', col_name_1)
						time_2 = re.findall('_([0-9]+\.?[0-9]?)_', col_name_2)
						if time_1 == time_2:
							col_1 = self.data[col_name_1]
							col_2 = self.data[col_name_2]
							correlation = col_1.corr(col_2)
							correlations.append(correlation)
				mean_corr = np.mean(correlations)
				if mean_corr > highest_mean_corr:
					highest_mean_corr = mean_corr
					cell_types = [pre_1, pre_2]
		return cell_types
		
	def most_static_genes(self, list_of_prefixes, number_of_genes):	
		"""
		Parameters: list of cell type prefixes, number of top genes to return
		Returns: list of genes with minimum average (across all genes) 
			standard deviation of expression across time points)
		"""
		gene_names = self.data['Gene Accession Number']
		self.data['std_dev'] = pd.Series(0, index=gene_names.index)
		# Clever to create a new col in the dataframe itself
		for row_index in range(0, len(self.data)):
			std_devs = []
			for prefix in list_of_prefixes:
				col_names = self.relevant_column_names(prefix)
				relevant_columns = self.data[col_names]
				relevant_exp_values = relevant_columns.ix[row_index]
				std_dev = np.std(relevant_exp_values)
				std_devs.append(std_dev)
			mean_std_dev = np.mean(std_devs)
			self.data.ix[row_index, 'std_dev'] = mean_std_dev
		sorted = self.data.sort('std_dev')
		sorted_gene_names = sorted["Gene Accession Number"]
		most_static_genes = sorted_gene_names[:number_of_genes]
		return most_static_genes
		
	def twofold_higher_genes_across_times(self, list_of_prefixes, time_1, time_2):
		"""
		Parameters: list of cell type prefixes, start time, end time
		Returns: list of genes whose expression at time 2 was at least twice as high as at time 1
		"""
		str_1 = str(time_1)
		str_2 = str(time_2)
		answer_genes = []
		col_names_1 = []
		col_names_2 = []
		for prefix in list_of_prefixes:
			col_name_1 = prefix + "_" + str_1 + "_hrs"
			col_name_2 = prefix + "_" + str_2 + "_hrs"
			col_names_1.append(col_name_1)
			col_names_2.append(col_name_2)
		for gene in range(0, len(self.data)):
			is_two_fold_different = 1
			for name in range(0, len(col_names_1)):
				name_1 = col_names_1[name]
				name_2 = col_names_2[name]
				at_time_1 = self.data[name_1][gene]
				at_time_2 = self.data[name_2][gene]
				fold_dif = (at_time_2 - at_time_1) / at_time_1
				if fold_dif < 2:
					is_two_fold_different = 0
			if is_two_fold_different == 1 :
				answer_genes.append(self.data["Gene Accession Number"][gene])
		return answer_genes
				
	def twofold_different_genes_across_two_cell_types(self, prefixes, time):
		"""
		Parameters: list containing two cell type prefixes, time
		Returns: list of genes whose expression at the specified time 
			exhibited at least a twofold difference between cell types
		"""
		answer_genes = []
		prefix_1 = prefixes[0]
		prefix_2 = prefixes[1]
		col_name_1 = prefix_1 + "_" + str(time) + "_hrs"
		col_name_2 = prefix_2 + "_" + str(time) + "_hrs"
		for gene in range(0, len(self.data)):
			in_cell_type_1 = self.data[col_name_1][gene]
			in_cell_type_2 = self.data[col_name_2][gene]
			fold_dif = (in_cell_type_2 - in_cell_type_1) / in_cell_type_1
			if math.fabs(fold_dif) > 2:
				answer_genes.append(self.data["Gene Accession Number"][gene])	
		return answer_genes
	

if __name__ == '__main__':
	expression_file = str(sys.argv[1])
	print "Bootcamp homework 1"
	ea = ExpressionAnalyzer()
	ea.import_file(expression_file)
	
	# Problem A
	print "A. \n\t" , ea.count_genes(), " genes"
	
	all_prefixes = ["HL60", "U937", "NB4", "Jurkat"]
	
	#Problem B
	tp_by_cell_type = ea.most_correlated_time_points(all_prefixes)
	print "B. "
	for cell_type in tp_by_cell_type:
		tp_pair = cell_type[1]
		print "\t", cell_type[0], ": ", tp_pair[0], " and ", tp_pair[1]
	
	#Problem C
	print "C. \n\t", ea.most_similar_cell_types(all_prefixes)
	
	#Problem D #(takes by far the longest)
	print "D. "
	i = 1
	for index, value in (ea.most_static_genes(all_prefixes, 10)).iteritems():
		print "\t", i, ": ", value
		i += 1
	
	#Problem E
	print "E. "
	twofold_genes_by_time = ea.twofold_higher_genes_across_times(all_prefixes, 0, 24)	
	for col in twofold_genes_by_time:
		print col
	
	#Problem F
	two_prefixes = ["HL60", "U937"]
	print "F. "
	twofold_genes_by_cell_type = ea.twofold_different_genes_across_two_cell_types(two_prefixes, 0)
	for col in twofold_genes_by_cell_type:
		print col
	
	#Problem G #could do directly through DAVID web client
	print "G. "
	
	
	