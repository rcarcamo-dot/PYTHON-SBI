#Import modules

import sequence_class 
import sequence_data
import re
import sys
import os
import argparse
import gzip
import random

#Define command-line arguments parser

parser = argparse.ArgumentParser(description="""This program """)

parser.add_argument('-i', '--input',
					dest = "infile",
					action = "store",
					default = None,
					help = "Input File")

parser.add_argument('-o', '--output-file',
					dest = "outfile",
					action = "store",
					default = None,
					help = "Output File")

parser.add_argument('-v', '--verbose',
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "Print log in stderr")


options = parser.parse_args()

## Function definitions


#Get the correspondent files of input

def get_input_file(input):
	""" Handling with different kind of input """
	fasp = re.compile('.pdb$|.ent$|.pdb.gz$|.ent.gz$')
	path = input

	if path == None:
		path = os.getcwd()

	if os.path.isdir(path):
		fasta_files = [f for f in os.listdir(path) if fasp.search(f) is not None]
		os.chdir(path)
	else:
		fasta_files = [input]

	return fasta_files

#Open different kind of output files

def get_output_file(output):
	""" Handling with different kind of output: None, normal file or gunzip file """

	if output == None:
		outfile = sys.stdout
	else:
		if ".gz" in output:
			outfile = gzip.open(output, "wt")
		else:
			outfile = open(output, "w")

	return outfile

def get_progress_log():
	"""
	Return a message in the standard error with the progression log of the program 
	"""
	pass

#Function to print data into outfile

def print_results(prot_list, outfile):
	"""
	Printing the results to selected outfile
	"""
	pass
#Execution Control from the original file:

if __name__ == "__main__":

	fasta_files = get_input_file(options.infile)
	outfile = get_output_file(options.outfile)

	info_prots = get_protein_len_mw(fasta_files)
	prots = info_prots[0]
	n_files = info_prots[1]
	n_sequences = info_prots[2]
	n_proteins = info_prots[3]

	if options.number:
		prots = get_n_proteins(prots, options.number)
		n_random = len(prots)
	else:
		n_random = 0

	if options.verbose:
		get_progres_log(n_files, fasta_files, n_sequences, n_proteins, n_random, options.outfile)

	print_results(prots, outfile)