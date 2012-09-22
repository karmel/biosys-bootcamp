from pandas import DataFrame

class GeneNumber(object):

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

	# convert a list of genes into a set to remove the duplicates, and get the length
	def countGene(self):
		count = len(set(self.dataframe['Gene Accession Number']))
		print("The number of genes is " + str(count))
		return

lst = GeneNumber([])
lst.import_file("data_set_HL60_U937_NB4_Jurkat.txt")
lst.countGene()
