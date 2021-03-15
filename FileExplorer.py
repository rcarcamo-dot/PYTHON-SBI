import Bio.PDB
import gzip
import io


""" Essentially what is happening here is that we are receiving a 
file along with what that file is going to be called and an indicator of 
whether or not it is zipped. Then this files is going to be sent to a PDBParser
which finds the chains located within the file and returns those chain objects. 
the chain objects contain information on the structure of the chain and can be used to 
find the sequence which in turn can be used to align. With this object we can also superimppose later on"""

class File(object):
	""" A class that takes as an input the path to a file 
	and returns various relevant info about it """

	def __init__(self, FilePath, name, Zip):
		self.FilePath = FilePath
		self.name = name
		self.Zip = Zip

	def get_FilePath(self):
		return self.FilePath

	def get_name(self):
		return self.name

class pdb(File):
	"""  A pdb child of the File contains 
	info on the chains found within the file """

	def get_chain(self):
		""" Returns a list of all of the chain onjects found within the provided file
		This is done with a PDBParser that can read both pdb and pdb.gzip files"""
		chainlist = []
		parser = Bio.PDB.PDBParser()
		if self.Zip == "FALSE" :
			PDBstructure = parser.get_structure(self.get_name(), self.get_FilePath())
			for model in PDBstructure:
				for chain in model:
					chainlist.append(chain)
			return chainlist

		else:
			with gzip.open(self.FilePath, 'rb') as finput:
				finput = io.StringIO(finput.read().decode("utf-8"))
				PDBstructure = parser.get_structure(self.get_name(), finput)
				for model in PDBstructure:
					for chain in model:
						chainlist.append(chain)
				return chainlist

