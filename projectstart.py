import sys
import argparse
import gzip
import os
import FileExplorer
import FunctionHodgepodge as FH
import functions

""" First we grab all of the user input as specified below"""

parser = argparse.ArgumentParser(description = "This program is the beginning of the Python/SBI project")
required = parser.add_argument_group('required arguments')

required.add_argument("-i","--input-directory",
	dest = "input",
	action = "store",
	default = None,
	required = True,
	help = "Input directory containing pdb files")

parser.add_argument("-s","--stoichiometry",
	dest = "stechiometry",
	action = "store",
	default = None,
	help = "input of stoichiometry in case a protein is a homodimer")

required.add_argument("-o","--output-directory",
	dest = "output",
	action = "store",
	default = None,
	help = "Input directory in which output files should be stored")

parser.add_argument("-f","--force",
	dest = "force",
	action = "store_true",
	default = False,
	help = "Input FASTA formatted file")

parser.add_argument("-v","--verbose",
	dest = "verbose",
	action = "store_true",
	default = False,
	help = "Print log in stderr")


""" now we save all of the arguments to options. We call them with optons"""
options = parser.parse_args()

""" here I am simply checing whether or not the input directory is an actual folder or if something went wrong.
If the input is not a directory I kill the program 
Carles has this as a seperate function I believe which we can implement"""

try:
	files = os.listdir(options.input)
except:
	print("You did not provide a valid input directory!")
	exit(0)

try:
	destination = os.listdir(options.output)
except:
	print("You did not provide a valid output directory!")
	exit(0)



Allinteractions = [] # This variable stores the name of the files within the input directory

####
###  The files storage variable above may have to be changed to account for other types of files contained besides just pfb (interactions)
####

for file in files:
	""" Here I am going to get each file that was provided in the directory 
	I willl simultaneously check if that file needs to be read with gzip or not"""
	i = 0
	split = file.split("_")  
	chainlist = [] # chain names are stored in a list that I can use to check if they appear in the file
	for line in split:
		line = line.replace(".pdb","")
		if ".gz" in line:
			line = line.replace(".gz","")
			i = 1
		if line != split[0]:
			chainlist.append(line.strip())


	""" Now I am going to go through to make sure that the chains indicated in the file header are actually found within the file
	First I look for all of the chain names and save this information"""
	filelocation = options.input + "/" + file
	try:
		if i == 1:
			with gzip.open(filelocation, "r") as fd:
				chainsfound = set()
				for row in fd:
					for chain in chainlist:
						if bytes(chain,encoding='utf8') in row:
							chainsfound.add(chain)
		else:					 
			with open(filelocation, "r") as fd:
				chainsfound = set()
				for row in fd:
					for chain in chainlist:
						if chain in row:
							chainsfound.add(chain)

		#print (chainsfound)
	except IsADirectoryError:
		print (f"The directory provided contains a directory!\n{filelocation}\n")

	""" Now I will check to make sure that the chain names in the title of the file match the chain names within the file itself """
	chainsgiven = set()
	for element in chainlist:
		chainsgiven.add(element)
	if len(chainlist) != len(chainsfound):
		print(chainsgiven.difference(chainsfound) , " was not found in the file")
		exit(0)


	Allinteractions.append(file)

### Now we are going to get chain information for each of the files in the directory 
AllChains = {} # We are going to save all of the chains here
for file in Allinteractions:
	filepath = options.input + "/" + file
	Zip = "FALSE"
	if ".gz" in file:
		Zip = "TRUE"
	pdb = FileExplorer.pdb(filepath, file, Zip )
	AllChains[file] = pdb.get_chain()
""" So after that last step we now have a list of chain objects from which we can retrieve information 
This list is stored as a value inside a dictionary that has the filename from which they came from as the key"""


""" Now I will try to align the chains across file pairs to try and find matches
In order to do this I will first have to get the sequences of each chain object so that I can align them."""

""" Here i am going to get the sequences for each chain and save them in the dictionary that I created"""
for file in AllChains:
	i=0
	for chain in AllChains[file]:
		if functions.get_molecule_type(chain) == "Protein":	
			AllChains[file][i] = functions.get_backbone_atoms_protein(chain)
		else:
			AllChains[file][i] =  functions.get_backbone_atoms_nucleicacids(chain)
		i += 1




#print (AllChains)



	#print(chainsgiven.difference(chainsfound))
	#print(len(chainlist))




