#!/usr/bin/env python3
"""
program to optimise positions of H atoms in a molecule
as read from an xyz file 
"""
import numpy as np
from ase import Atoms

from copy import deepcopy

import ase_tools
import sphere_point_picking

def get_boltzmann(dE, RT=0.2):
    """
    get Boltzmann term exp(- dE / RT)
    """
    exponent = -1 * dE / RT
    return np.exp(exponent)

def get_random_vector():
    """
    get vector pointing in random direction and with random 
    length
    """
    # get length sampled from normal distribution with mu=0.5
    # sigma=0.05
    r = np.random.normal(0.5, 0.05)
    # get corresponding coords 
    coords = sphere_point_picking.sample_sphere_cartesian(r)
    return coords

def update_coords(*coordinates):
    """
    add random vector to all coordinates 
    """
    for coord in coordinates: 
        coord += get_random_vector()
    return coordinates

def update_atoms(atoms):
    # don't modify old atoms obj in memory
    nuclear_charges = deepcopy(atoms.numbers)
    coordinates = deepcopy(atoms.positions)

    # H atoms to be modified
    H_indices = np.where(nuclear_charges == 1)[0]
    
    # change H coordinates
    # randomly choose which to modify 
    index = np.random.choice(H_indices, size=1)
    coordinates[index] = update_coords(coordinates[index])
    new_atoms = ase_tools.init_atoms_obj(nuclear_charges, coordinates)
    return new_atoms

def optimise_geometry(atoms, maxiter=1000, RT=0.25, deldump=False):
    """
    using PM7 in MOPAC
    WARNING: MOPAC will dump a bunch of files: set deldump=True
    to delete these

    Parameters
    ----------
    atoms: ase.Atoms object 
           Initial geometry to optimise 
    maxiter: int 
             number of iterations of MC optimisation
    RT: float 
        temperature for the Boltzmann distribution
    deldump: Bool 
             whether or not to delete MOPAC's output files

    Returns
    -------
    atoms_new: ase.Atoms object
               Final geometry
    E_new: float 
           Energy of final geometry
    """
    # initialise variables
    atoms_old = atoms
    E_old = ase_tools.get_energy(atoms_old)
    print("original energy", E_old)

    # counter
    accept = 0

    # begin simulation 
    for i in range(maxiter):
        # translation move 
        atoms_new = update_atoms(atoms_old)

        # calc energy of new geometry 
        E_new = ase_tools.get_energy(atoms_new)
        # sometimes doesn't return energy (?)
        if E_new:
            dE = E_new - E_old

            # Metropolis-Hastings acceptance criterion
            boltzmann = get_boltzmann(dE, RT)
            if boltzmann >= np.random.uniform(): 
                # accept new geometry
                atoms_old = atoms_new 
                E_old = E_new
                accept += 1

           # else:
                # leave as is
            #    print("reject new geometry")

    acceptance_ratio = (accept / maxiter) * 100
    print("acceptance ratio", acceptance_ratio, "over", maxiter, "iterations")

    if deldump:
        ase_tools.del_outfiles()

    return atoms_new, E_new


if __name__ == "__main__":
    trial = "../data/dsgdb9nsd_000004/path_0_frag_0_fragment.xyz"
    atoms = ase_tools.init_atoms_from_file(trial)
    atoms_new, E_new = optimise_geometry(atoms, maxiter=1000, deldump=True)
    print("final energy", E_new)
    # save file
#    ase_tools.write_atoms_to_file(atoms_new, "opt.xyz")










