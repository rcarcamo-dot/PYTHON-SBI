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

A python package to construct macrocomplex given protein or nucleic acid pairwise interactions.

## Tutorial
### Input Files

The input files must be <b>pairs of interacting chains</b> (.pdb), which has to be located into a directory decided by the user. This directory can also be in <i>tar.gz</i> format.

### Python Modules

A module is a file containing Python definitions and statements. It is important to consider that definitions from a module can be imported into other modules or into the <i>main</i> module.

* <b>sbi_project.py</b>: <i>main</i> module or program created to reconstruct a macrocomplex given a set of interacting pairs (prot-prot, prot-RNA). It also considers the possibility to rebuild a macrocomplex if the input is a macrocomplex (.pdb). 

  Aditionally, this module contains the ArgumentParser object, created using argparse module, which is used for command-line     options, arguments and sub-commands. The ArgumentParser object will hold all the information necessary to parse the command   line into Python data types.
  
* <b>functions.py</b>: this module is composed by a set of different <b>functions</b> to solve biological and technical problems during the analysis. Thus, it is imported into the other modules in order to use the defined functions.  
  
  






