from Bio.PDB import Superimposer

class Superimposer (object):
    """
    """
    def __init__ (self, object1, object2):
        self.object1 = object1
        self.object2 = object2
        self.super_imposer = PDB.Superimposer()
    

    def imposer (self):
        
        if molecule_type == "Protein":
            atoms1 = list(self.object1.get_backbone_atoms_protein)
        else:
            atoms2 = list (self.object2.get_backbone_atoms_nucleicacids )

        if len(atoms1) > len(atoms1):
            atoms1 = atoms1[:len(atoms2)]
        else:
            atoms2 = atoms2[:len(atoms1)]
        
        self.Superimposer.set_atoms(atoms1, atoms2)

        return self.Superimposer

    def calculate_RMSD (self)

        imposer = self.imposer
        rmsd = imposer.rms
        return rmsd

    print(super_imposer.rms)