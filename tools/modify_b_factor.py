from Bio.PDB.MMCIFParser import MMCIFParser


parser = MMCIFParser()

structure = parser.get_structure("struct_name", "file.cif")

for atom in structure.get_atoms():
    if atom.name == "CA":
        idx += 1
        atom.set_bfactor(value_to_map[idx])
    else:
        atom.set_bfactor(0.0)


from Bio.PDB.mmcifio import MMCIFIO
io=MMCIFIO()
io.set_structure(structure)
io.save("struct_with_modfied_bfactor.cif")
