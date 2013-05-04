from pandas import DataFrame

class ExpressionAnalyzer(object):

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

	def moreThan2Fold(self):
		# list of sample's name. More samples can be added very easily
		lst = ['HL60_0_hrs', 'HL60_24_hrs', 'U937_0_hrs', 'U937_24_hrs', 'NB4_0_hrs', 'NB4_24_hrs', 'Jurkat_0_hrs', 'Jurkat_24_hrs']
		count = 0
		samples = []
		fl = open('2FoldGenes.txt','w') 	# open a file for write
		#retrieve sample
		for name in lst:
			samples.append(self.dataframe[name])

		for i in range(0, len(samples[0])):
			j=0
			while(True):
				ratio = samples[j+1][i] / samples[j][i]
				j = j + 2
				if(ratio < 2 or j+1 > len(samples)):
					break
			if(ratio > 2):
				fl.write(self.dataframe.ix[i][0] + '\n')
				count += 1
		fl.write('\n' + str(count) + " Genes")
		print('File "2FoldGenes.txt" has been created!\n')
		return

lst = ExpressionAnalyzer([])
lst.import_file("data_set_HL60_U937_NB4_Jurkat.txt")
lst.moreThan2Fold()
