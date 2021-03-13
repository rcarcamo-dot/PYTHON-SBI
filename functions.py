
import sys
import os
import re
from Bio.PDB import PDBParser



class IncorrectInput(NameError):
	def __init__(self,input):
		self.input = input
	def __str__(self):
		return "The directory name %s is not a directory " %(self.input)

def obtain_pdb_files (directory):
	""" This function obtains the ".pdb" files from a specified directory and returns a list of the aforementioned files. """

	directory = os.listdir(directory) # Saves all the files in a variable.
	list_of_pdb_files = []
	for file in directory:
		if file.endswith(".pdb"):
			list_of_pdb_files.append(file) # Saves the pdb files in a list.
	return list_of_pdb_files


def obtain_structure(pdb_name,pdb_file):
    """ This function parses a PDB file to obtain a structure object. """

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_file[0:-4], pdb_file)
    return structure

def get_backbone_atoms_protein (chain):
	""" This function returns a list of proteins' backbone atoms (CA) from a given chain. """

	ca_atoms = []
	for residue in chain:
		if residue.get_id()[0] == " " and residue.has_id("CA"):
			ca_atoms.append(residue['CA'])
	return ca_atoms

def get_backbone_atoms_nucleicacids (chain):
	""" This function returns a list of nucleic acids' backbone atoms (C4') from a given chain. """

	c4_atoms = []
	for residue in chain:
		if residue.get_id()[0] == " " and residue.has_id("C4\'"):
			c4_atoms.append(residue["C4\'"])
	return c4_atoms

def get_molecule_type (chain):
	""" This function identifies the molecule type (RNA, DNA, PROTEIN) of a given chain and returns it as a string. """
	# Creates a list for each DNA and RNA chain with all possible letter for eachone.
	RNA = ['A','U','C','G','I']
	DNA = ['DA','DT','DC','DG','DI']
	molecule_type = ""

	for residue in chain:
		residue_name = residue.get_resname().strip()
		break
	if residue_name in RNA:
		molecule_type = "RNA"
	elif residue_name in DNA:
		molecule_type = "DNA"
	else:
		molecule_type = "Protein"
	return molecule_type


def get_best_core(list_of_pdb_files):
	""" This function selects the file containing the chain with larger number of interactions from the list of input files to act as the core for the superimposition strategy, and returns this file back. """

	recount = {}
	## Creates a dictionary with the different chains in input_directory files and count the times they appear ##
	for file in list_of_pdb_files:
		if file[-5] in recount:
			recount[file[-5]] += 1
		if file [-6] in recount:
			recount[file[-6]] += 1
		if file[-5] not in recount:
			recount.setdefault(file[-5], 1)
		if file[-6] not in recount:
			recount.setdefault(file[-6], 1)
	recount_sorted = {k: v for k, v in sorted(recount.items(), key=lambda item: item[1], reverse=True)} # Sorts the dictionary by its values, which represent the most frequent chain in pdb files from input_directory.
	chains = list(recount_sorted.keys())
	core_chain = chains[0] # Select the most frequent chain.
	for file in list_of_pdb_files:
		search_pattern = re.compile(core_chain) # Saves the chain id as a search pattern.
		m = search_pattern.search(file, 5, 7) # Search the most frequent chain through all pdb files
		if m:
			core_chain_file = file # Select the first pdb file contaning the search pattern as core chain file.
			break

	return core_chain_file