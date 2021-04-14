import sys
import argparse
import gzip
import os
import FileExplorer
import functions
import superimposer
import Bio.PDB
import random
import copy
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.PDB import Superimposer

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
	type = int,
	dest = "stoichiometry",
	action = "store",
	default = None,
	help = "input of stoichiometry in case a protein is a homo-mer (number of identical chains in complex)")

required.add_argument("-o","--output-directory",
	dest = "output",
	action = "store",
	default = None,
	help = "Input directory in which output files should be stored")

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



##################                  ##################
##################  Initial Setup   ##################
##################                  ##################

Allinteractions = [] # This variable stores the name of the files within the input directory
File_chain_pair = [] # This variable contains tuples of all of the chains and files that they are located in
####
###  The files storage variable above may have to be changed to account for other types of files contained besides just pdb (interactions)
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


# AllChains is now a list of tuples containing the chains in the second position and the file wherin that chai can be found in the first

if options.verbose:
	print("checking input and output locations...\n\n")

newchain_id = 1
newchain_id_old_id = {}

### Now we are going to get chain information for each of the files in the directory 
AllChains = {} # We are going to save all of the chains here
for file in Allinteractions:
	AllChains[file] = {}
	filepath = options.input + "/" + file
	Zip = "FALSE"
	if ".gz" in file:
		Zip = "TRUE"
	# I send the information to the pdb extractor.
	pdb = FileExplorer.pdb(filepath, file, Zip )
	contained_chain_names = []
	for element in file.split(".pdb")[0].split("_"):
		if file.split(".pdb")[0].split("_").index(element) != 0:
			contained_chain_names.append(element)
	

	""" What is happening here is the creation of dictionary that will allow us to access a chain 
	object with the filename and the name of the chain. For this I am setting up a new unique id for
	each chain that will allow me to use it as key in a dictionary and then having a translation dictionary
	that can let me know which original chain that belonged to"""
	""" I was originally going to create a subclass of the chain class that contained a hashing function
	so that I could just use the chain as a key, but since I though of this after a bunch of coding this solution 
	was easier than going back to rewrite the class that I used. """

	for chain in pdb.get_chain():
		onetuple = (file,chain)
		File_chain_pair.append(onetuple)
		old_id = chain.id
		AllChains[file][old_id] = {}
		chain.id = newchain_id
		newchain_id_old_id[newchain_id] = (file, old_id)
		newchain_id += 1
		AllChains[file][old_id][chain.id] = ""
		#exit(0) 


####
####
####
""" this following loop prints the translation dictionary """

#for newchain_id in newchain_id_old_id:
#	print (newchain_id, "->    ", newchain_id_old_id[newchain_id])

""" Here I can see exacty what the ids are"""
#for file in AllChains:
#	for old_id in AllChains[file]:	
#		print(AllChains[file][old_id])
#		for new_id in AllChains[file][old_id]:
#			print(newchain_id_old_id[new_id])
####
####
####

""" So after that last step we now have a list of chain objects from which we can retrieve information 
This list is stored as a value inside a dictionary that has the filename from which they came from as the key"""



##################                             ##################
##################   Information Processing    ##################
##################                             ##################


""" Now I will try to align the chains across file pairs to try and find matches
In order to do this I will first have to get the sequences of each chain object so that I can align them."""

if options.verbose:
	print("getting Ca atoms for superpositioning...\n\n")

# new_id_to_chain_dict = {}  <<< incase i need this (superimposer does it)

""" Backbone atoms are retrieved later now. Here I am just saving the chain. Incase this changes in future iterations of 
the program I am keeping the structure to save the backbone atoms commented"""
for file in AllChains:
	i=0
	for chain in AllChains[file]:
		for chain_object_id in AllChains[file][chain]:
			for filecheck, chaincheck in File_chain_pair:
				if filecheck == file and chaincheck.id == chain_object_id:
					#new_id_to_chain_dict[chaincheck.id] = chaincheck <<< incase i need this (superimposer does it)
					if functions.get_molecule_type(chaincheck) == "Protein":
						# file is the first key (file of origin for the chain)
						# chain is second key (original chain id)
						# chain_object_id is third key (new chain id)
						# I am now adding the chain as the final value
						AllChains[file][chain][chain_object_id] = chaincheck
						#AllChains[file][chain][chain_object_id] = functions.get_backbone_atoms_protein(chaincheck)  <<< incase i need this (superimposer does it)
					else:
						AllChains[file][chain][chain_object_id] =  chaincheck
						#AllChains[file][chain][chain_object_id] =  functions.get_backbone_atoms_nucleicacids(chaincheck) <<< incase i need this (superimposer does it)
					i += 1	


""" in order to avoid calculating the sequences multiple times
I am just going place the sequences and their identifiers in a list of tuples to be used
in the Alignment"""

if options.verbose:
	print("Saving sequences...\n\n")

forAlignmentlist = []
for file in AllChains:
	i = 0
	for chain in AllChains[file]:
		for new_id in AllChains[file][chain]:
			addtuple = (file + "_" + chain, functions.get_sequence(AllChains[file][chain][new_id]) )
			forAlignmentlist.append(addtuple)
			i += 1




#############                            #############
#############   Dealing with Homo-mers   #############
#############                            #############


""" If the protein is a Homo-mer we will only have a single file to deal with 
in that case the user has to provide us with an integer of the amount of repeating chains
of that homomer in the end structure."""

###
### The method of growing the model here is the same as the original which will occur 
### Further on (if no stoichiometry was provided). The process is explained there
###

# checking if the user gave us any stoichiometry input
if options.stoichiometry != None:
	# variables set up to control for the amount of chains that we are going to add to the model
	iterations = options.stoichiometry - 2
	number_of_homos_added = 0
	model_indicator = 1
	working_model = Bio.PDB.Model.Model(model_indicator)
	Homomer_chain_list = []
	model_chain_id_contained = {}
	for file in AllChains:
		for original_id in AllChains[file]:
			for new_id in AllChains[file][original_id]:
				working_model.add(AllChains[file][original_id][new_id])
				Homomer_chain_list.append(AllChains[file][original_id][new_id])
				model_chain_id_contained[AllChains[file][original_id][new_id].id] = 1

	if options.verbose:
		print("the first two chains of the homo-mer were added to the model...\n\n")
	while number_of_homos_added < iterations:
		for fixed_chain in working_model.get_chains():
				moving_chain = Homomer_chain_list[1]
				growth_chain  = copy.copy(Homomer_chain_list[0])
				# I add an id to the new growth chain that will be generated at random	
				new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				while new_growth_chain.id in model_chain_id_contained:
					new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				for residue in growth_chain:
					new_growth_chain.add(residue.copy())
				working_model.add(new_growth_chain)
				number_of_homos_added += 1

				fixed_atoms = functions.get_the_atoms(fixed_chain)
				moving_atoms = functions.get_the_atoms(moving_chain)
				growth_atoms = functions.get_the_atoms(growth_chain)

				if len(fixed_atoms) != len(moving_atoms):
					worked_chains = functions.get_worked_chains(fixed_chain, moving_chain)
					worked_fixed_chain = worked_chains[0]
					worked_moving_chain = worked_chains[1]

					fixed_atoms = functions.get_the_atoms(worked_fixed_chain)
					moving_atoms  = functions.get_the_atoms(worked_moving_chain)

				superpositioning = Bio.PDB.Superimposer()
				superpositioning.set_atoms(fixed_atoms,moving_atoms)
				superpositioning.apply(new_growth_chain.get_atoms())

				flag = True
				if not functions.ensure_no_clashing(working_model, new_growth_chain):
					working_model.detach_child(new_growth_chain.id)
					number_of_homos_added = number_of_homos_added - 1
					flag = False
				if flag:
					model_chain_id_contained[new_growth_chain.id] = 1
					if options.verbose:
						print("The Protein welcomes a new chain into the fold...\n\n")



	manualOveride = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k',
              'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F',
              'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '!',
              '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', ':', ';', '<', '=', '>', '?', '@',
              '[', ']', '^', '_', '`', '{', '|', '}', '~']

	for x in working_model:
		x.id = manualOveride.pop(0)
		#print(x.id)
	#print(model_chain_id_contained)
	for file in files:
		name = file.split("_")[0]
	model_filename = options.output +"/"+ name + ".pdb"
	io = Bio.PDB.PDBIO()
	io.set_structure(working_model)
	io.save(model_filename)

	if options.verbose:
		print("###\n### The PDB file of the protein has been saved at: ", model_filename, "\n###\n\n")



	exit(0)


if options.verbose:
	print("Aligning chains to search for matches...\n\n")
""" Here I am going to perform the alignments. i am saving the results of the 
alignment to a file named Resuts_of_Alignments.txt that will be located in the user
provided outputs directory. In order to identify the iddentical chains I take the alignment
score and devide it by the maximum sequence length of the two aligned sequences.This normalizes the score.
Now I can apply the user provided threshold to find identical sequences."""
writeresult = open(options.output + "/Results_of_Alignments.txt", "w")
Filesdone = {}
chainpair_dict = {}
writeresult.write("id1\tid2\tscore\n")
for chain in forAlignmentlist:
	if options.verbose:
		print(f"Aligning chain {chain[0]}...\n")
	Filesdone[chain[0]] = ""
	for secondchain in forAlignmentlist: 
		if secondchain[0] not in Filesdone and chain[0][:-2] != secondchain[0][:-2]:
			correctinglength = max(len(chain[1]),len(secondchain[1]))
			Align_result = pairwise2.align.globalxx(chain[1],secondchain[1])
			if not len(Align_result) == 0:
				writeresult.write(chain[0]+"\t"+secondchain[0]+"\t"+str(Align_result[0][2]/correctinglength)+"\n")
				if options.threshold <=  Align_result[0][2]/correctinglength:
					if chain[0] in chainpair_dict:
						chainpair_dict[chain[0]].append(secondchain[0])
					else:
						chainpair_dict[chain[0]] = [secondchain[0]]
				#print(Align_result[0][2])

if options.verbose:
	print("###\n### Alignments are saved in outputs directory as \"Results_of_Alignments.txt\"\n###\n\n")
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



##################                        ##################
##################   Model Construction   ##################
##################                        ##################


if options.verbose:
	print("Building up the first chains of the model...\n\n")
""" Here I am constructing the starting structure of the model with the intitial reference pair of chains
that I calculated in the previous step. The information for the chains contained in the model is saved is 
stored in the model_chain_id_contained dictionary (by the orginal chain id) this is to know what further matches
I can add to the model. Only interactions that contain a chain that is found in the model are interesting."""
# a variable that saves what chain is in the model (need this to check for matches)
model_chain_id_contained = {}
# a modelindicator serves as the id, I could take user input to act as id as well. Not necessary though
model_indicator = 1
# The working model starts here (called working model because it isnt the completed one yet)
working_model = Bio.PDB.Model.Model(model_indicator)
# getting the file of the initial chains
reference_file = initialref.split(".pdb_")[0] + ".pdb"
# adding the two chains to the model. We now have the beginnings of a model
working_model.add(list(AllChains[reference_file][initialref.split(".pdb_")[0].split("_")[1]].values())[0])
model_chain_id_contained[list(AllChains[reference_file][initialref.split(".pdb_")[0].split("_")[1]].values())[0].id] = 1
working_model.add(list(AllChains[reference_file][initialref.split(".pdb_")[0].split("_")[2]].values())[0])
model_chain_id_contained[list(AllChains[reference_file][initialref.split(".pdb_")[0].split("_")[2]].values())[0].id] = 1

#exit(0)
""" Here I am constructing a dictionary that saves the information about the matches that were found.
all of the chains involved in the relationship are saved here with the role being reflected in the 
position in the tuple. I also save the static alternate chain incase I need it (its part of the 
failsafe that I implement in the model constructon)the information is saved in a tuple acting as a key in the dictonary """
if options.verbose:
	print("Setting up matches ...\n\n")
chain_match_list = {}
for x in chainpair_dict:
	static_file = x.split(".pdb_")[0] + ".pdb"
	static_chain = list(AllChains[static_file][x.split(".pdb_")[1]].values())[0]
	for char in x.split(".pdb_")[0].split("_"):
		if x.split(".pdb_")[0].split("_").index(char) != 0 and char != x.split(".pdb_")[1]:
			static_alt_chain = list(AllChains[static_file][char].values())[0]

	for string in chainpair_dict[x]:
		chain_name = string.split(".pdb_")[1]
		file = string.split(".pdb_")[0] + ".pdb"
		for char in string.split(".pdb_")[0].split("_"):
			if string.split(".pdb_")[0].split("_").index(char) != 0 and char != chain_name:
				alt_chain = list(AllChains[file][char].values())[0]
		chain = list(AllChains[file][chain_name].values())[0]
		chain_match_list[(static_chain,chain,alt_chain,static_alt_chain)] = 1

""" Here I am going through the matches and checking to find what I can line up. I set a while loop 
that iterates as long as there are still fixed chain references to go by. I break the while loop as 
soon as there have been no new chain additions to the model because if nothing new was added then the next 
run through will give me no new results. If however a single chain is added to the model I must continue looping through """
# here I store all of the matches that have already been made (this is just for trouble shooting)
if options.verbose:
	print("Beginning Model Construction...\n\n")
tuples_used = set()
# variable storing the number of additions per round
number_of_addition = 1
while number_of_addition != 0:
	number_of_addition = 0	
	# I go through each chain in the model (when I add a new one I go through it in the next while loop)
	for fixed_chain in working_model.get_chains():
		# now I match that chain up with the moving chains
		for match_tuple in chain_match_list:
			if match_tuple[0] == fixed_chain:
				moving_chain = match_tuple[1]
				growth_chain  = copy.copy(match_tuple[2])
				# I add an id to the new growth chain that will be generated at random	
				new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				while new_growth_chain.id in model_chain_id_contained:
					new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				for residue in growth_chain:
					new_growth_chain.add(residue.copy())
				working_model.add(new_growth_chain)
				number_of_addition += 1

				fixed_atoms = functions.get_the_atoms(fixed_chain)
				moving_atoms = functions.get_the_atoms(moving_chain)
				growth_atoms = functions.get_the_atoms(growth_chain)

				if len(fixed_atoms) != len(moving_atoms):
					worked_chains = functions.get_worked_chains(fixed_chain, moving_chain)
					worked_fixed_chain = worked_chains[0]
					worked_moving_chain = worked_chains[1]

					fixed_atoms = functions.get_the_atoms(worked_fixed_chain)
					moving_atoms  = functions.get_the_atoms(worked_moving_chain)

				superpositioning = Bio.PDB.Superimposer()
				superpositioning.set_atoms(fixed_atoms,moving_atoms)
				superpositioning.apply(new_growth_chain.get_atoms())

				flag = True
				if not functions.ensure_no_clashing(working_model, new_growth_chain):
					working_model.detach_child(new_growth_chain.id)
					number_of_addition = number_of_addition - 1
					flag = False
				if flag:
					model_chain_id_contained[match_tuple[2].id] = 1
					model_chain_id_contained[new_growth_chain.id] = 1
					if options.verbose:
						print("The Protein welcomes a new chain into the fold...\n\n")
				tuples_used.add(match_tuple)
				
		""" Because I arbitrarily decided that the first chain in the alignments would 
		be considered the 'static' chain I have to now check to make sure that if I used
		the second chain isn't included in the model for some reason 
		THis is just a failsafe and should not be necessary based on how I set up the alignment
		(in the sense that I only compared going forward and didn't align back. I may be able
		to remove this repetition"""


		for match_tuple in chain_match_list:
			if match_tuple[1] == fixed_chain:
				moving_chain = match_tuple[0]
				growth_chain  = copy.copy(match_tuple[3])
				new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				while new_growth_chain.id in model_chain_id_contained:
					new_growth_chain = Bio.PDB.Chain.Chain(random.randint(0,100000))
				for residue in growth_chain:
					new_growth_chain.add(residue.copy())
				working_model.add(new_growth_chain)
				number_of_addition += 1

				fixed_atoms = functions.get_the_atoms(fixed_chain)
				moving_atoms = functions.get_the_atoms(moving_chain)
				growth_atoms = functions.get_the_atoms(growth_chain)

				if len(fixed_atoms) != len(moving_atoms):
					worked_chains = functions.get_worked_chains(fixed_chain, moving_chain)
					worked_fixed_chain = worked_chains[0]
					worked_moving_chain = worked_chains[1]

					fixed_atoms = functions.get_the_atoms(worked_fixed_chain)
					moving_atoms  = functions.get_the_atoms(worked_moving_chain)

				superpositioning = Bio.PDB.Superimposer()
				superpositioning.set_atoms(fixed_atoms,moving_atoms)
				superpositioning.apply(new_growth_chain.get_atoms())

				flag = True
				if not functions.ensure_no_clashing(working_model, new_growth_chain):
					working_model.detach_child(new_growth_chain.id)
					number_of_addition = number_of_addition - 1
					flag = False
				if flag:
					model_chain_id_contained[match_tuple[3].id] = 1
					model_chain_id_contained[new_growth_chain.id] = 1
					if options.verbose:
						print("The Protein welcomes a new chain into the fold...\n\n")
				tuples_used.add(match_tuple)





if options.verbose:
	print("Construction has finished...\n\n")




""" I had an issue using integers as chain ids and found out that they have to be characters so I set up this overide to fix it """
### THIS IS THE ISSSUE the ids have to be characters!!!!!!!!!
manualOveride = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k',
              'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D', 'E', 'F',
              'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '!',
              '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', ':', ';', '<', '=', '>', '?', '@',
              '[', ']', '^', '_', '`', '{', '|', '}', '~']

for x in working_model:
	x.id = manualOveride.pop(0)
	#print(x.id)
#print(model_chain_id_contained)

name = reference_file.split("_")[0]
model_filename = options.output + "/" + name + ".pdb"
io = Bio.PDB.PDBIO()
io.set_structure(working_model)
io.save(model_filename)

if options.verbose:
	print("###\n### The PDB file of the protein has been saved at: ", model_filename, "\n###\n\n")



exit(0)





########
########
########
########
########  THE END 
########
########
########
########
########










