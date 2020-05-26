#!/usr/bin/env python3
"""
test case for the optimiser
"""
import ase_tools
import optHs

def ethene_distorted(savefile=False):
    """
    this should be quite easy to optimise as 
    I only slightly distorted one of the C and 
    H positions in the molecule
    """
    file = "../data/dsgdb9nsd_000004/sample.xyz"
    atoms = ase_tools.init_atoms_from_file(file)
    atoms_new, E_new = optHs.optimise_geometry(atoms, maxiter=1000,
                                                deldump=True)
    print("final energy", E_new)

    if savefile:
        ase_tools.write_atoms_to_file(atoms_new, "opt.xyz")

if __name__ == "__main__":
    ethene_distorted()
