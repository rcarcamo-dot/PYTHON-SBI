
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


def process_stoichiometry:

    stoich = options.stoichiometry
    stoich_dict = {}
    i = 0
    j = 0
    while i < len(stoich):
        if stoich[i].isalpha():
            stoich_dict[stoich[i]] = ""
            j += 1
            while stoich[j].isdigit():
                stoich_dict[stoich[i]] += str(stoich[j])
                if j < (len(stoich) - 1):
                    j += 1
                else:
                    break
        i += 1
        j = i

    stoich_dict = {x: int(y) for x, y in stoich_dict.items()}
    k = sum(list(stoich_dict.values()))

    if options.verbose:
        sys.stderr.write("Stoichiometry: %s\n" % options.stoichiometry)
	
def check_stoichiometry:
    stoich_chains = set(list(similar_chains.values()))

    if len(stoich_chains) < len(stoich_dict):
        raise ValueError("Unmatching stoichiometry: Provided stoichiometry contains %d different chains, but "
                         "the input PDBs only have %d unique chains." % (len(stoich_dict), len(stoich_set)))
