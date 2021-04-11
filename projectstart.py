import sys
import argparse
import gzip
import os
import FileExplorer
import FunctionHodgepodge as FH
import functions
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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
	help = "######")

parser.add_argument("-v","--verbose",
	dest = "verbose",
	action = "store_true",
	default = False,
	help = "Print log in stderr")

parser.add_argument("-t", "--threshold",
	type = float, 
	dest = "threshold",
	action = "store",
	default = float(0.9),
	help = "Alignment score threshold for sequences to be considered indentical")


""" now we save all of the arguments to options. We call them with optons"""
options = parser.parse_args()

""" here I am simply checing whether or not the input directory is an actual folder or if something went wrong.
If the input is not a directory I kill the program 
Carles has this as a seperate function I believe which we can implement"""

if options.verbose:
	print("checking input and output locations...\n\n")

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

if options.verbose:
	print("Checking format of input files and gathering file and chain names...\n\n")

Allinteractions = [] # This variable stores the name of the files within the input directory
AllChains = [] # This variable contains tuples of all of the chains and files that they are located in
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
	for chain in chainlist:
		onetuple = (file,chain)
		AllChains.append(onetuple)

# AllChains is now a list of tuples containing the chains in the second position and the file wherin that chai can be found in the first

if options.verbose:
	print("checking input and output locations...\n\n")

### Now we are going to get chain information for each of the files in the directory 
AllChains = {} # We are going to save all of the chains here
for file in Allinteractions:
	AllChains[file] = {}
	filepath = options.input + "/" + file
	Zip = "FALSE"
	if ".gz" in file:
		Zip = "TRUE"
	pdb = FileExplorer.pdb(filepath, file, Zip )
	for chain in pdb.get_chain():
		AllChains[file][chain] = {}
""" So after that last step we now have a list of chain objects from which we can retrieve information 
This list is stored as a value inside a dictionary that has the filename from which they came from as the key"""


""" Now I will try to align the chains across file pairs to try and find matches
In order to do this I will first have to get the sequences of each chain object so that I can align them."""

if options.verbose:
	print("getting Ca atoms for superpositioning...\n\n")

""" Here i am going to get the Ca atoms for eventual superposition for each chain and save them in the dictionary that I created"""
for file in AllChains:
	i=0
	for chain in AllChains[file]:
		if functions.get_molecule_type(chain) == "Protein":	
			AllChains[file][chain] = functions.get_backbone_atoms_protein(chain)
		else:
			AllChains[file][chain] =  functions.get_backbone_atoms_nucleicacids(chain)
		i += 1


""" in order to avoid calculating the sequences multiple times
I am just going place the sequences and their identifiers in a list of tuples
or maybe a dictionary"""

if options.verbose:
	print("Saving sequences...\n\n")

forAlignmentlist = []
for file in AllChains:
	i = 0
	for chain in AllChains[file]:
		addtuple = (file + "_" + chain.id, functions.get_sequence(chain) )
		forAlignmentlist.append(addtuple)
		i += 1

if options.verbose:
	print("Aligning chains to search for matches...\n\n")
""" Here I am going to perform the alignments. i am saving the results of the 
alignment to a file named Resuts_of_Alignments.txt that will be located in the user
provided outputs directory. In order to identify the iddentical chains I take the alignment
score and devide it by the maximum sequence length of the two aligned sequences.
It is probably a better idea to filter out the identical sequences now rather than going throug
this file later to get the matches. I would use a user provided threshold to calculate this."""
writeresult = open(options.output + "/Results_of_Alignments.txt", "w")
Filesdone = {}
chainpair_dict = {}
writeresult.write("id1\tid2\tscore\n")
for chain in forAlignmentlist:
	Filesdone[chain[0]] = ""
	for secondchain in forAlignmentlist: 
		if secondchain[0] not in Filesdone and chain[0][:-2] != secondchain[0][:-2]:
			correctinglength = max(len(chain[1]),len(secondchain[1]))
			Align_result = pairwise2.align.globalxx(chain[1],secondchain[1])
			if options.threshold <=  Align_result[0][2]/correctinglength:
				writeresult.write(chain[0]+"\t"+secondchain[0]+"\t"+str(Align_result[0][2]/correctinglength)+"\n")
				if chain[0] in chainpair_dict:
					chainpair_dict[chain[0]].append(secondchain[0])
				else:

					chainpair_dict[chain[0]] = [secondchain[0]]
				#print(Align_result[0][2])

if options.verbose:
	print("###\n### Alignments are saved in outputs folder as \"Results_of_Alignments.txt\"\n###\n\n")

""" So the matching Chains have been identified, and are stored in a file within the outputs directory.
The file where the information can be found is called "Results_of_Alignments.txt". The format of the file
is tab seperated. In the first two columns are the identifiers of the two chain pairs. The way how the chains are
identified are by "their filename" + "_" + "the chain name". The information can be extracted for superpositioning."""

if options.verbose:
	print("Creating a file to place the models in if it doesn't exist yet...\n\n")

created_models = options.output  + "/created_models"
if not os.path.exists(created_models):
	os.makedirs(created_models)

if options.verbose:
	print("Checking for the best starting reference chains...\n\n")
##
## Turn this into a function
##

""" I am simply checking which of the chains has the most matches associated with it 
It is not necessary, but we need to start somewhere and I figure this is a good place to do so """
initialref = ""
mostmatch = 0
for key in chainpair_dict:
	if len(chainpair_dict[key]) > mostmatch:
		initialref = key
		mostmatch = len(chainpair_dict[key])

if options.verbose:
	print("###\n### The best starting reference chain is", initialref, " with ",  mostmatch ,"matches\n###\n\n")


##
## that is function
##

##
## this could also be a function
##
for match in chainpair_dict[initialref]:
	reference_file = initialref.split(".pdb_")[0] + ".pdb"
	reference_chain = initialref.split(".pdb_")[1]
	match_file = match.split(".pdb_")[0] + ".pdb"
	match_chain = match.split(".pdb_")[1]
	for element in match.split(".pdb_")[0].split("_"):
		if match.split(".pdb_")[0].split("_").index(element) != 0 and element != reference_chain:
			growth_chain = element
			print(growth_chain)

	print (initialref, " \t", match)
	print (reference_file, " ", reference_chain,"\t", match_file, " ", match_chain,"\n\n")
##
## that is function
##







