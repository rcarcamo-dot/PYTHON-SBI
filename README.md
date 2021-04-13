# PYTHON-SBI
## Constructing Macrocomplexes

*Patrick Gohl, Oumout Egkemen Moustafa, Roberto Carcamo Calvo*

## **Table of Contents**

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

## Prerequisites

* <b>Python</b> 3.0 Release

https://www.python.org/download/releases/3.0/

* <b>Modeller</b> 10.1 Release

https://salilab.org/modeller/download_installation.html

* <b>Chimera</b> 1.15 Release

https://www.cgl.ucsf.edu/chimera/download.html

## Theoretical Background

Most of the proteins comprise of a variety of polypeptide chains that make up multi-subunit complexes. Such form is known as the quaternary structure of a protein, which plays a significant role in determining its function, activity and localization. Additionally, these protein subunits may be interacting with nucleic acids (DNA) and ribonucleic acids(RNA) as well. However, comprehending how individual subunits form into bigger complexes is not an easy work. In this project, we created a program to create multi-subunit complexes out of its individual pairwise interactions. 

## Tutorial
### Input Files

The input files must be pairs of interacting chains (.pdb), which has to be located into a directory decided by the user. This directory can also be in <i>tar.gz</i> format.

### Python Modules

A module is a file containing Python definitions and statements. It is important to consider that definitions from a module can be imported into other modules or into the <i>main</i> module.

* <b>projectstart.py</b>: <i>main</i> module created to reconstruct a macrocomplex given a set of interacting pairs. It also considers the possibility to rebuild a macrocomplex if the input is a macrocomplex (.pdb). 

  Furthermore, this module inludes the ArgumentParser object, created using argparse module, which is used for command-line options, arguments and sub-commands.
  
* <b>functions.py</b>: this module is composed by a set of different <b>functions</b> to solve biological and technical problems during the analysis. Thus, it is imported into the other modules in order to use the defined functions.  

* <b>superimposer.py</b>: this module includes codes, functions and libraries imported in order to superimpose given a set of interacting pairs.

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
    
* <b> -s, --stoichiometry</b>

  This option will allow user to input of stoichiometry in case a protein is homodimer. 

  <i>(default: None)</i>  
  
* <b> -t, --threshold</b>

  This option will allow user to determine threshold for sequences to be considered identical. 

  <i>(default: 0.9)</i>  


If the <b>default options</b> are set, these are the following outputs:

<b>files</b>: 

 - Result_of_Alignments.txt: file genereated from the pairwise comparisons between all the inputs.



