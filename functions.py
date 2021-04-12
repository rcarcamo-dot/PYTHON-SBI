import sys
import os
import re
from Bio.SeqUtils import IUPACData as threetoone 
from Bio.PDB import PDBParser
from Bio.PDB import Chain
import Bio
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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


""" get sequence is supposed to return the full sequence of the chain"""


def get_sequence (chain):
	sequence = ""
	for residue in chain:
		if residue.get_id()[0] == " " and residue.has_id("CA"):
			addition = threetoone.protein_letters_3to1[residue.get_resname().capitalize()]
			sequence = sequence + addition

	return sequence


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


def ensure_no_clashing(model, growth_chain):
	""" ensures that the new chain isnt interfering (clashing) with any of the existing 
	chains that in out model. Returns True if it isnt and false if it does"""
	clashing = set()
	model_atom_list = Bio.PDB.Selection.unfold_entities(model, "A")
	Neighbors = Bio.PDB.NeighborSearch(model_atom_list)

	for atoms in growth_chain.get_atoms():
		location = atoms.get_coord()
		matches = Neighbors.search(location, 1.2, level = "C")
		for result in matches:
			if result != growth_chain:
				clashing.add(result)
	if len(clashing) == 0:
		return True
	else:
		return False

def prep_superimpose(fixed_chain, moving_chain):
	fixed_chain = fixed_chain
	moving_chain = moving_chain
	fixed_chain_sequence = get_sequence(fixed_chain) 
	moving_chain_sequence = get_sequence(moving_chain)

	alignments = pairwise2.align.globalxx(fixed_chain_sequence, moving_chain_sequence)

	fixed_chain_sequence = alignments[0][0]
	moving_chain_sequence = alignments[0][1]

	fixedchain_pattern = []
	movingchain_pattern = []

	for res1, res2 in zip(fixed_chain_sequence, moving_chain_sequence):
		if res1 == res2:
			fixedchain_pattern.append(1)
			movingchain_pattern.append(1)
		elif res1 == '-':
			movingchain_pattern.append(0)
		elif res2 == '-':
			fixedchain_pattern.append(0)
	chains_pattern = (fixedchain_pattern, movingchain_pattern)
	return chains_pattern

def get_worked_chains(fixed_chain, moving_chain):
	new_fixed_chain = Bio.PDB.Chain.Chain('X')
	new_moving_chain = Bio.PDB.Chain.Chain('Y')
	fixed_chain = fixed_chain
	moving_chain = moving_chain
	for residue, pattern in zip(fixed_chain.get_residues(), prep_superimpose(fixed_chain, moving_chain)[0]):
		if pattern == 1:
			new_fixed_chain.add(residue.copy())
	for residue, pattern in zip(moving_chain.get_residues(), prep_superimpose(fixed_chain, moving_chain)[1]):
		if pattern == 1:
			new_moving_chain.add(residue.copy())
	return (new_fixed_chain,new_moving_chain)

def get_the_atoms(chain):
	if get_molecule_type(chain) == "Protein":
		atoms = list(get_backbone_atoms_protein(chain))
	else:
		atoms = list (get_backbone_atoms_nucleicacids(chain))
	return atoms






