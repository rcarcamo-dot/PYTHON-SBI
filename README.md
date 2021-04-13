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

A program to construct macrocomplex given interacting pairs (protein-protein, protein-nucleic acid, protein - RNA).

## Prerequisites

* <b>Python</b> 3.0 Release

https://www.python.org/download/releases/3.0/

* <b>Modeller</b> 10.1 Release

https://salilab.org/modeller/download_installation.html

It is also recommended to have a program for interactive visualization and analysis of molecular structures and related data such as <b>Chimera</b> .

* <b>Chimera</b> 1.15 Release

https://www.cgl.ucsf.edu/chimera/download.html

## Theoretical Background

Most of the proteins comprise of a variety of polypeptide chains that make up multi-subunit complexes. Such form is known as the quaternary structure of a protein, which plays a significant role in determining its function, activity and localization. Additionally, these protein subunits may be interacting with nucleic acids (DNA and RNA) as well. However, comprehending how individual subunits form into bigger complexes is not an easy work. In this project, we created a program to create multi-subunit complexes out of its individual pairwise interactions.

In the process of modelling complexes, one of the fundemental steps is to find out which chains shows up in more than one input pair. To evaluate it, we made decision for performing a pairwise sequence allignment between all the chains from the whole set of input pairs. Since we thought that working with normalized values would be better for modeling, the values were divided by the length of the longest chain in the alignment. We considered chains with normalized values 0.9 and higher as equal.

After finding out chains of interest, the latter step is to place the chains of our complexes into three dimensional space. Since superimposition is one of the most used method for this purpose, we used the module Superimposer of Bio.PDB. In this module the input chains are into three categories : fixed, moving and growth chain. The chain produces the maximum number of matches in the file is called as "fixed" and when we perform superimposition, it will not be moving. The chains will be superimposed of the fixed chain are called moving chains. One crucial thing in this point is that fixed and moving chains need to have the equal atom length. In case there are some gaps in the alignment, we have to remove those residues to make sure that list of atoms completely lines up. Superimpose applies a rotation and translation matrix onto every atom of the moving chain to put into the same position as the fixed chain. This rotation and translation movements will be kept in the rotran matrix of the Superimposer object. Applying the rotran matrix to the growth chain which does not match up with the first building block results in a change of its coordinates and now it functions like fixed chain. Insertion of growth chain may lead to occurance  of steric clashes. For this reason when a new chain added to model, we ensure no clashing by checking that no atom in the growth chain after applying the rotran matrix is occupying the same spaces as the rest of model. In case steric clash is detected, we remove the chain and we do not use it. This process will be iterated as many times as needed. As the models grows, there will be more and more chains to check for matches. We finish the build when no chain has been added.

## Tutorial
### Input Files

The input files must be pairs of interacting chains (.pdb), which has to be located into a directory decided by the user. This directory can also be in <i>tar.gz</i> format.

### Python Modules

A module is a file containing Python definitions and statements. It is important to consider that definitions from a module can be imported into other modules or into the <i>main</i> module.

* <b>projectstart.py</b>: <i>main</i> module created to reconstruct a macrocomplex given a set of interacting pairs. It also considers the possibility to rebuild a macrocomplex if the input is a macrocomplex (.pdb). Furthermore, this module inludes the ArgumentParser object, created using argparse module, which is used for command-line options, arguments and sub-commands.  
  
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
    
* <b> -s, --stoichiometry</b>

  This option will allow user to input of stoichiometry in case a protein is homodimer. 

  <i>(default: None)</i>  
  
* <b> -t, --threshold</b>

  This option will allow user to determine threshold for sequences to be considered identical. 

  <i>(default: 0.9)</i>  


If the <b>default options</b> are set, these are the following outputs:

<b>files</b>: 

 - Result_of_Alignments.txt: file genereated from the pairwise comparisons between all the inputs.



