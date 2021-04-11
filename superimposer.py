from Bio.PDB import Superimposer
import functions
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
        else:
            atoms1 = list (functions.get_backbone_atoms_nucleicacids(self.object1))

        if functions.get_molecule_type(self.object2) == "Protein":
            atoms2 = list(functions.get_backbone_atoms_protein(self.object2))
        else:
            atoms2 = list (functions.get_backbone_atoms_nucleicacids(self.object2))

        if functions.get_molecule_type(self.growth_chain) == "Protein":
            atoms3 = list(functions.get_backbone_atoms_protein(self.growth_chain))
        else:
            atoms3 = list (functions.get_backbone_atoms_nucleicacids(self.growth_chain))


        if len(atoms1) > len(atoms2):
            atoms1 = atoms1[:len(atoms2)]
        else:
            atoms2 = atoms2[:len(atoms1)]
        
        Superimposer.set_atoms(self, atoms1, atoms2)

        Superimposer.apply(self, self.growth_chain.get_atoms())

        return self.growth_chain

    def get_growth_chain(self):
        return self.growth_chain