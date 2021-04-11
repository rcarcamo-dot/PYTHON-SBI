sp1 = []
sp2 = []

for model in object1:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                    sp1.append(atom)

for model in object2:
    for chain in model:
        if chain.get_id() == sp_chain:
            for residue in chain:
                for atom in residue:
                    sp2.append(atom)


super_imposer = PDB.Superimposer()
super_imposer.set_atoms(sp1, sp2)
super_imposer.apply(object2.get_atoms())


print(super_imposer.rms)


io = PDB.PDBIO()
io.set_structure(object2)
io.save("s2.pdb")
io.set_structure(object1)
io.save("s1.pdb")