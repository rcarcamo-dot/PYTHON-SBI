# PYTHON-SBI
## Constructing Macrocomplexes

*Patrick Gohl, Oumout Egkemen Moustafa, Roberto Carcamo Calvo*

## **TABLE OF CONTENTS**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [PYTHON-SBI](#python-sbi)
  - [Table of Contents](#Table-of-Contents)
  - [Description](#Description)
- [Prerequisites](#Prerequisites)
- [Installation](#Installation)
- [Theoretical Background](#Theoretical-Background)
- [Tutorial](#Tutorial)
  - [Input Files](#Input-Files)
  - [Python Modules](#Python-Modules)
  - [Output Files](#Output-Files)
  - [Analysis of Reconstructed Macromolecules](#Analysis-of-Reconstructed-Macromolecules)
- [Limitations](#Limitations) 
  
<!-- /TOC -->


## Description

A python package to construct macrocomplex given interacting pairs (protein-protein, protein-nucleic acid, protein - RNA).

## Tutorial
### Input Files

The input files must be <b>pairs of interacting chains</b> (.pdb), which has to be located into a directory decided by the user. This directory can also be in <i>tar.gz</i> format.

### Python Modules

A module is a file containing Python definitions and statements. It is important to consider that definitions from a module can be imported into other modules or into the <i>main</i> module.

* <b>projectstart.py</b>: <i>main</i> module created to reconstruct a macrocomplex given a set of interacting pairs. It also considers the possibility to rebuild a macrocomplex if the input is a macrocomplex (.pdb). 

  Aditionally, this module contains the ArgumentParser object, created using argparse module, which is used for command-line options, arguments and sub-commands. The ArgumentParser object will hold all the information necessary to parse the command line into Python data types.
  
* <b>functions.py</b>: this module is composed by a set of different <b>functions</b> to solve biological and technical problems during the analysis. Thus, it is imported into the other modules in order to use the defined functions.  

### Output Files

It is necessary to consider that the output files will depend on the command-line options and arguments the user determines while performing the analysis. All of them will be stored in the output directory selected by the user (-o or --output-directory).

The following <b>arguments</b> are established:

* <b> -i, --input-directory  </b>
  
  It must be a directory provided by the user which contains the inputs pairs (compressed format is also available <i>.tar.gz</i>).

   <i>(default: None)</i>
 

* <b> -o, --output-directory</b>
 
  This is a directory which will be created during the program, structured in other subdirectories. 
  
  <i>(default:None)</i>
  
  
* <b> -v, --verbose</b>

  This option will allow the user to completely follow the progam. 

  <i>(default: False)</i>
  
  
* <b> -f, --force</b>

  This option will allow user to input FASTA formatted file. 

  <i>(default: False)</i>
  
  
* <b> -s, --stoichiometry</b>

  This option will allow user to input of stoichiometry in case a protein is homodimer. 

  <i>(default: None)</i>  
  

* <b> -t, --threshold</b>

  This option will allow user to determine threshold for sequences to be considered identical. 

  <i>(default: 0.9)</i>  


If the <b>default options</b> are set, these are the following outputs:

<b>files</b>: 

 - Result_of_Alignments.txt: file genereated from the pairwiise comparisons between all the inputs.



