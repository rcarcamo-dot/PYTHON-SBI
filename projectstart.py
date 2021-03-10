import sys
import argparse
import gzip
import os
import FileExplorer

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



options = parser.parse_args()

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

AllChains = [] # We are going to save all of the chains here
for file in Allinteractions:
	filepath = options.input + "/" + file
	Zip = "FALSE"
	if ".gz" in file:
		Zip = "TRUE"
	pdb = FileExplorer.pdb(filepath, file, Zip )
	AllChains.append(pdb.get_chain())


print (AllChains)



	#print(chainsgiven.difference(chainsfound))
	#print(len(chainlist))




