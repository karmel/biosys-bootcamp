from pandas import DataFrame
from numpy import array
from heapq import nsmallest

class GetGenes(object):

	def __init__(self, data):
		self.dataframe = DataFrame(data)

	# read a text file and return a data frame. Records should be separated by TAB
	# There should not be duplicate column names
	def import_file(self, filename):
		# this function use to convert string to float
		def convert(x):
			try:
				x = float(x)
			except ValueError:
				pass
			return(x)

		table = []
		for line in open(filename):
			if(line.strip()):	# If not empty line
				line = line.rstrip('\n').split('\t')
				line = list(map(convert, line))
				table.append(line)
		self.dataframe = DataFrame(table[1:],columns=table[0])
		return

	def houseKeepingGenes(self, geneNum):
		# compute the CV of data
		std = array(self.dataframe.std(axis = 1))
		mean = array(self.dataframe.mean(axis = 1))
		CV = std/mean
		CV = list(map(abs, CV))		# convert to positive number

		# get the fist N minimum value
		mins = nsmallest(geneNum, CV)
		print("The GOOD genes are:\n")
		for item in mins:
			print(self.dataframe.ix[CV.index(item)][0])
		return

lst = GetGenes([])
lst.import_file("data_set_HL60_U937_NB4_Jurkat.txt")
lst.houseKeepingGenes(10)
