from pandas import DataFrame

class DiffExpression(object):

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

	def diffGenes(self, s1, s2):
		fl = open('differentiated_expression.txt', 'w')
		for i in range(0, len(self.dataframe[s1])):
			# Make sure to import division from __future__ if you do float division, just to be safe
			ratio = self.dataframe[s1][i] / self.dataframe[s2][i]
			if(ratio >= 2 or ratio <= 0.5):
				fl.write(self.dataframe.ix[i][1] + '\n')
		print('Your file "differentiated_expression.txt" has been written!\n')
		return

lst = DiffExpression([])
lst.import_file("data_set_HL60_U937_NB4_Jurkat.txt")
lst.diffGenes('HL60_0_hrs', 'U937_0_hrs')
