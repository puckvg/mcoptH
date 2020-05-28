#!/usr/bin/env python3
"""
This file calculates the energy of a molecule with 
particular atomtypes and coordinates using MOPAC
within ASE
"""

from ase import Atoms 
from ase import io
from ase.calculators.mopac import MOPAC
from ase.optimize import LBFGS

import os

def init_atoms_obj(nuclear_charges, coordinates):
    atoms = Atoms(numbers=nuclear_charges, positions=coordinates)
    return atoms 

def init_atoms_from_file(filename):
    atoms = io.read(filename)
    return atoms 

def write_atoms_to_file(atoms, filename, format="xyz", append=False):
    io.write(filename=filename, images=atoms, format=format, append=append)

def get_energy(atoms, calculator=MOPAC, task="1SCF", method="PM7"):
    # initialise calculator
    calc = calculator(label="mopac", method=method, task=task)
    # assign to atoms obj
    atoms.calc = calc 

    # calc energy - sometimes for some reason this isn't possible, 
    # so we return None
    try:
        e = atoms.get_potential_energy()
    except:
        return None

    return e

def opt(atoms, calculator=MOPAC, task="1SCF GRADIENTS", method="PM7", opt_method=LBFGS):
    atoms.calc = calculator(label="mopac", task=task, method=method)
    dyn = opt_method(atoms)
    dyn.run()

    del_outfiles()

    e = atoms.get_potential_energy()
    return atoms, e

def del_outfiles():
    """
    delete MOPAC dump - assuming the files were labelled with mopac
    """
    os.remove("mopac.arc")
    os.remove("mopac.mop")
    os.remove("mopac.out")




