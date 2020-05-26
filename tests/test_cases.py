#!/usr/bin/env python3
"""
test case for the optimiser
"""
import ase_tools
import optHs

def good_guess():
    """
    good initial guess of ethene
    """
    file = "good_initial_guess.xyz"
    atoms = ase_tools.init_atoms_from_file(file)
    atoms_new, E_new = optHs.optimise_geometry(atoms, maxiter=1000,
                                               deldump=True)
    ase_tools.write_atoms_to_file(atoms_new, "opt_ethene.xyz")
    print("final energy", E_new)

def bad_guess():
    """
    this is a terrible initial guess of methane
    """
    file = "bad_initial_guess.xyz"
    atoms = ase_tools.init_atoms_from_file(file)
    atoms_new, E_new = optHs.optimise_geometry(atoms, maxiter=1000,
                                                deldump=True)
    ase_tools.write_atoms_to_file(atoms_new, "opt_methane.xyz")
    print("final energy", E_new)

if __name__ == "__main__":
    good_guess()
    bad_guess()
