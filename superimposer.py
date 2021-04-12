from Bio.PDB import Superimposer
import functions
from Bio.PDB import NeighborSearch
import numpy


class Superposer (object):
    """
    """
    def __init__ (self, object1, object2, growth_chain):
        self.object1 = object1
        self.object2 = object2
        self.growth_chain = growth_chain

    def imposer (self):

        
        if functions.get_molecule_type(self.object1) == "Protein":
            atoms1 = list(functions.get_backbone_atoms_protein(self.object1))
            chain1 = (object1, atoms1) #save chains in a list, check if it's better to store them in a hash
        else:
            atoms1 = list (functions.get_backbone_atoms_nucleicacids(self.object1))
            chain1 = (object1, atoms1)

        if functions.get_molecule_type(self.object2) == "Protein":
            atoms2 = list(functions.get_backbone_atoms_protein(self.object2))
            chain2 = (object1, atoms2)
        else:
            atoms2 = list (functions.get_backbone_atoms_nucleicacids(self.object2))
            chain2 = (object2, atoms2)

        if functions.get_molecule_type(self.growth_chain) == "Protein":
            atoms3 = list(functions.get_backbone_atoms_protein(self.growth_chain))
            chain3 = (object3, atoms3)
        else:
            atoms3 = list (functions.get_backbone_atoms_nucleicacids(self.growth_chain))
            chain3 = (object3, atoms3)



        if len(object1[1]) > len(object2[1]):
            object1[1] = object1[1][:len(obeject2[1])]
        elif len(object2[1]) > len(object1[1]):
            object2[1] = object2[1][:len(obeject1[1])]
        
        return object2
        


    def get_homodimer(imposer):
        """
        """
        homodimers[]
        if object1 == object2 :
            assert numpy.abs(superimposition.rms) < 0.000000001
            assert numpy.max(numpy.abs(superimposition.rotran[1])) < 0.00000001
            assert numpy.max(numpy.abs(superimposition.rotran[0]) - numpy.identity(3)) < 0.00000001
            homodimers[(chain2, atoms2)]
        else:
            pass
           
        print("RMS(first model, model %i) = %0.2f" % (object2.id, superimposition.rms))

        return homodimer


    def complex_builder (self)
        """
        """
        superimposition = Superimposer() 
			superimposition.set_atoms(object1[1], object2[1])
			RMSD = superimposition.rms 
			if RMSD > RMSD_threshold: 
				continue
            else:
                pass
			superimpositions[(object1.get_id(),object2.get_id())] = superimposition #check if we already have a get_id function
	        if superimposed_chains is True: 
	            superimpositions = sorted(superimpositions.items(), key = lambda x: x[1].rms) 
	    return (superimpositions, best_RMSD, superimposed_chains)