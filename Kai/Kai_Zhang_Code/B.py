from pandas import DataFrame

class MaxCorrelation(object):

	# An interesting approach to make separate modules for each part of the question.
	# If you do so, though, you should make a parent class so that you don't have to repeat
	# reused code, like the init function, in each file.
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

	# take a list of column names, compute the correlation by a pairwise way, and print the pair with maximum correlation
	def maxCorr(self, nameLst):
		n = len(nameLst)
		maxCo = 0
		sampleName = []
		for i in range(0, n):
			for j in range(i+1, n):
				corr = self.dataframe[nameLst[i]].corr(self.dataframe[nameLst[j]])
				if(abs(corr) > abs(maxCo)):
					maxCo = corr
					sampleName = [nameLst[i], nameLst[j]]
		print("Most highly correlated samples in ", nameLst,"are " + sampleName[0] + " AND " + sampleName[1] + ". Value: " + str(maxCo))
		return # Unlike R, in Python, you don't have to have return unless you actually want to return something.
		# The default is that the function returns None, so you don't need return here.

# Sample Name List
nameLst1 = ['HL60_0_hrs', 'HL60_0.5_hrs', 'HL60_4_hrs', 'HL60_24_hrs']
nameLst2 = ['U937_0_hrs', 'U937_0.5_hrs', 'U937_4_hrs', 'U937_24_hrs']
nameLst3 = ['NB4_0_hrs', 'NB4_5.5_hrs', 'NB4_24_hrs', 'NB4_48_hrs', 'NB4_72_hrs']
nameLst4 = ['Jurkat_0_hrs', 'Jurkat_0.5_hrs', 'Jurkat_4_hrs', 'Jurkat_24_hrs']

lst = MaxCorrelation([])
lst.import_file("data_set_HL60_U937_NB4_Jurkat.txt")
lst.maxCorr(nameLst1)
lst.maxCorr(nameLst2)
lst.maxCorr(nameLst3)
lst.maxCorr(nameLst4)
