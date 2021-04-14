# PYTHON-SBI
## Constructing Macrocomplexes

*Patrick James Gohl, Oumout Egkemen Moustafa, Roberto Carcamo Calvo*

## **Table of Contents**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [PYTHON-SBI](#python-sbi)
  - [Table of Contents](#Table-of-Contents)
- [Description](#Description)
- [Prerequisites](#Prerequisites)
- [Theoretical Background](#Theoretical-Background)
- [Tutorial](#Tutorial)
  - [Installation](#Installation)
  - [Input Files](#Input-Files)
  - [Python Modules](#Python-Modules)
  - [Command Line Arguments](#Command-Line-Arguments)
  - [Analysis of Reconstructed Macromolecules](#Analysis-of-Reconstructed-Macromolecules)
- [Limitations](#Limitations) 
  
<!-- /TOC -->


## Description

A program to construct macrocomplex given interacting pairs (protein-protein, protein-nucleic acid).

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

In the process of modelling complexes, one of the fundemental steps is to find out which chains shows up in more than one input pair. An assumption stated in the assignment of this project was that we would receive pdb files of all interacting chains of the protein that we are trying to reconstruct. So in order to reconstruct the complete protein from these files we simply had to connect the various pairs to eachother. To evaluate it, we made decision for performing a pairwise sequence allignment between all the chains from the whole set of input pairs. Since we thought that working with normalized values would be better for modeling, the values were divided by the length of the longest chain in the alignment. We considered chains with normalized values 0.9 and higher as equal. Once we know what chains are identical across all of the input files we can begin to reconstruct the whole complex one addition at a time.

After finding out chains of interest, the latter step is to place the chains of our complexes into three dimensional space. Since superimposition is one of the most used method for this purpose, we used the module Superimposer of Bio.PDB. In the process of superimposition we will have three chains of relevance. The fixed chain is the chain that will act as the reference. In our case it will always be the chain that is part of the model that we are adding to. The moving chain is the chain that produced a significant alignment with the fixed chain. THis means that it is the same chain. The growth chain as we are calling it is that chain which we know interacts with the moving chain, and it is the one that we want to ultimately add to the module. One crucial thing in this point is that fixed and moving chains need to have the equal atom length. In case there are some gaps in the alignment, we have to remove those residues to make sure that list of atoms completely lines up. Superimpose applies a rotation and translation matrix onto every atom of the moving chain to put into the same position as the fixed chain. This rotation and translation movements will be kept in the rotran matrix of the Superimposer object. Applying the rotran matrix to the growth chain which does not match up with the first building block results in a change of its coordinates. The new coordinates will place it into the proper interacting position in relation to the growing model that we are creating. Insertion of growth chain may lead to occurance  of steric clashes. For this reason when a new chain is added to model, we ensure no clashing by checking that no atom in the growth chain after applying the rotran matrix is occupying the same spaces as the rest of model. In case steric clash is detected, we remove the chain and we do not use it. This process will be iterated as many times as needed. As the models grows, there will be more and more chains to check for matches. We finish the build when no chain has been added.

Stoichiometry is another important aspect to be taken into account when we deal with homomers.Homomers present their own unique challenge to protein complex reconstruction. In this case the program will only be dealing with a single file that contains the information about one example of how the chains will interact with each other. Because the program would not be able to reconstruct the protein with this information alone it requires additional input in the way of stoichiometry. Once our program was able to handle this last hurdle we knew that it was ready for submission. In our program, the user has to provide us with an integer of the amount of repeating chains of that homomer in the end structure.

## Tutorial

### Installation

### Input Files

The input files must be pairs of interacting chains (.pdb), which has to be located into a directory decided by the user. This directory can also be in <i>tar.gz</i> format.

### Python Modules

A module is a file containing Python definitions and statements. It is important to consider that definitions from a module can be imported into other modules or into the <i>main</i> module.

* <b>projectstart.py</b>: <i>main</i> module created to reconstruct a macrocomplex given a set of interacting pairs. It also considers the possibility to rebuild a macrocomplex if the input is a macrocomplex (.pdb). Furthermore, this module inludes the ArgumentParser object, created using argparse module, which is used for command-line options, arguments and sub-commands.  
  
* <b>functions.py</b>: this module is composed by a set of different <b>functions</b> to solve biological and technical problems during the analysis. Thus, it is imported into the other modules in order to use the defined functions.  

* <b>FileExplorer.py</b>: This class is implemented to apply the PDB.Parser package in order to convert chains within a pdb file into chain objects. It is able to do this regardless of whether the files are gziped or not. 

### Command Line Arguments

It is necessary to consider that the output files will depend on the command-line options and arguments the user determines while performing the analysis. All of them will be stored in the output directory selected by the user (-o or --output-directory).

The following <b>arguments</b> are established:

* <b> -i, --input-directory  </b>
  
   The input directory should be an existing directory on you computer that contains pdb files of all the interacting chains of a protein that you whish to reconstruct. The command you would give would be  <b>(-i /PATH/to/the/directory)</b>. The name of the files should follow the format of given example : 1gzx_A_B.pdb


   <i>(default: None)</i>
 
* <b> -o, --output-directory</b>
 
   This is a directory which will be created during the program, structured in other subdirectories. The output directory is simply the path to the directory that you wish the output of the program to be saved to <b>(-o /PATH/to/output/folder)</b>. This must be an existing directory on your computer. 
  
  <i>(default:None)</i>
  
* <b> -v, --verbose</b>

   This is a command that tells the program to keep you updated as it is running the program. Just submit (-v) and the program will print to the terminal what exactly it is doing within the compiler. 

  <i>(default: False)</i>
    
* <b> -s, --stoichiometry</b>

   In case a protein is homodimer the stoichiometry command will allow you to tell the program how many instances of the chain occur in the protein that is to be reconstructed. For example if the protein is a Homo-4-mer, you will have to submit (-s 4) 

  <i>(default: None)</i>  
  
* <b> -t, --threshold</b>

  This option will allow user to set the threshold at which two sequences are considered equal after they have been aligned and had their score normalized. We recommend leaving the threshold at 0.9 as reducing it will create noise and cause the program to have to run through a lot of unnecessary calculations. However if you are expecting some of your sequences to be degraded in relation to others it may be a good idea to drop the threshold level. In this case you can submit (-t 0.8) for example. 

  <i>(default: 0.9)</i>  


### Analysis of Reconstructed Macromolecules

In order to analyze our program, we tested it against 3 different types of complexes. In order to analyze the quality of our program we decided that we would use the RMSD score. Since each of the protein that we worked with had a document 3D structure available on PDB we were able to simply Superimpose predicted and known structure in Chimera and calculate the RMSD. In order to superimpose, we used  matchmaker and chose to conduct chain pairing by means of the best aligning pair of chains between the known and predicted structure. We will go over the quality of the test structures again as their quality will also limit the quality of our analysis.

### 1GZX (the Simple Protein Complex)
The very first hurdle was being able to get the program to reconstruct a very simple complex that contained only amino acid sequences and was relatively short. As an example, we used the protein provided in the Python Aula Global “1gzx”. The problems presented by this protein were the ones shared throughout. These are mentioned in the Theory section. This protein acted as a sort of introduction for us on the problem of protein reconstruction, since it was so basic. 1gzx is an oxy state hemoglobin. The method by which the structure on pdb was obtained was X-Ray diffraction. The Resolution of that structure was 2.1 Angstrom. This is a pretty good resolution for a protein complex to have. Even for X-Ray diffraction it is around the average good result that is obtained. The quality of our reconstruction was perfect, achieving an RMSD of 0. The chimera reply log can be viewed in our submission under 1gzx_RMSD. We are also submitting the folder we used to reconstruct this complex (1gzx). 

Example command for this protein:
```python
python3 projectstart.py -i /Users/patrick/PYTHON/PythonProject/examples/1gzx -o /Users/patrick/PYTHON/outputs -t 0.9 -v
```

<p align="center">
  <img src="/1gzx.png" width="450"/>
</p>

### 5FJ8 (the Complex Protein)
The next step was to ensure that our program would work on a much more complicated protein structure. This one would contain significantly more chains as well as be bound to nucleic acid sequences instead of just other amino acid chains. To test our ability to reconstruct such a protein we selected the protein 5fj8 also supplied on the Python Aula Global web page. 5fj8 is a RNA Polymerase 3 elongation complex. On pdb it was obtained by means of Electron Microscopy at a resolution of 3.9 Angstroms. This is the lowest quality resolution that we will be using, but it isn’t that bad for EM. Again, after reconstructing the protein we superimposed the predicted and known structure in chimera. A disclaimer, because this protein complex is so large MatchMaker takes a long time to conduct the aligning to find the best chain pair. If you want to reconstruct the superimposing in a quicker fashion you could select a specific chain to use. The RMSD score obtained from the superimposition was perfect, at 0. We have included the chimera reply log under 5fj8_RMSD. We are also submitting the folder we used to reconstruct this complex (5fj8).

Example command for this protein:
```python
python3 projectstart.py -i /Users/patrick/PYTHON/PythonProject/examples/5fj8 -o /Users/patrick/PYTHON/outputs -t 0.9 -v 
```

<p align="center">
  <img src="/5fj8.png" width="450"/>
</p>

### 6C70 (The Homo-mer)
Finally, we needed to be able to reconstruct a homomeric protein. In order to tackle this challenge we had to create our own files to use since none of this type were provided in the Python Aula Global. We simply searched for a known homomer on PDB. We selected 6c70 (Orco) which is a Homo-4-mer. It is a Membrane protein and acts as an odorant receptor. The structure on PDB was obtained by means of Cryo - Electron Microscopy at a resolution of 3.5 Angstroms. For a protein obtained by Cryo-EM this is a good resolution. It is not however as good a resolution as one would get by x-ray crystallography which can yield between 1.5 – 2.5 Angstroms in the high resolution end. Nonetheless, for the purposes of testing our structure this is adequate. We selected this because it counterbalances the higher resolution proteins used previously as well as applying a different method to obtain the structure than the other ones. We downloaded the pdb file from PDB and worked it up in gedit text edit. All we had to do was remove all chains but A and B which were interacting. This file we saved and submitted to our program. After reconstructing the protein we Superimposed predicted and known structures in Chimera as before and again obtained a perfect RMSD score of 0. We are submitting the reply log as 6c70_RMSD as well as the file used to reconstruct this protein. 

Example command for this protein:
```python
python3 projectstart.py -i /Users/patrick/PYTHON/PythonProject/examples/6c70 -o /Users/patrick/PYTHON/outputs -t 0.9 -v -s 4
```

<p align="center">
  <img src="/6c70.png" width="450"/>
</p>

