#!/usr/bin/env python3
"""
program to optimise positions of H atoms in a molecule
as read from an xyz file 

Intended usage:
./optHs.py -xyz_in -n_steps -xyz_out
where xyz_in is the file to optimise, n_steps is the number of MC steps and xyz_out
is the outfile of the optimised structure 
"""
from copy import deepcopy
import argparse as ap

import numpy as np
from ase import Atoms

import ase_tools
import sphere_point_picking

def parse_args():
    parser = ap.ArgumentParser(description="read file and set options for optimisation")
    parser.add_argument("xyz_in", metavar="xyz_in", type=str, help="infile")
    parser.add_argument("n_steps", metavar="n", type=int, help="n optimisation steps")
    parser.add_argument("xyz_out", metavar="xyz_out", type=str, help="outfile")
    parser.add_argument("save_gif", metavar="gif", type=str, help="save gif yes/y or no/n")

    args = parser.parse_args()
    return args.xyz_in, args.n_steps, args.xyz_out, args.save_gif

def get_boltzmann(dE, RT=0.1):
    """
    get Boltzmann term exp(- dE / RT)

    Parameters 
    ----------
    dE: float 
        change in energy in kcal/mol
    RT: float 
        temperature in K multiplied by gas constant

    Returns 
    -------
    boltzmann: float 
               exp(- dE / RT)
    """
    exponent = -1 * dE / RT
    boltzmann = np.exp(exponent)
    return boltzmann

def get_random_vector(mu=0., sigma=0.2):
    """
    get vector pointing in random direction and with random 
    length

    Parameters 
    ----------
    mu: float 
        centre of normal distribution
    sigma: float 
           stdev of normal distribution

    Returns 
    -------
    coords: np.array 
            coordinates of a vector generated by choosing a random
            point on the surface of a sphere, and scaling by 
    """
    r = np.random.normal(mu, sigma)
    # get corresponding coords 
    coords = sphere_point_picking.sample_sphere_cartesian(r)
    return coords

def update_coords(*coordinates, mu=0., sigma=0.2):
    """
    add random vector to all coordinates 

    Parameters 
    ----------
    coordinates: list, np.array, list of lists, etc 
                 coordinates (vectors) to modify
    mu: float 
        centre of normal distribution
    sigma: float 
           stdev of normal distribution

    Returns
    -------
    coordinates: list, np.array, list of lists, etc 
                 matches input type 
                 coordinates (vectors) with random vectors added
    """
    for coord in coordinates: 
        coord += get_random_vector(mu, sigma)
    return coordinates

def update_atoms(atoms, mu=0., sigma=0.2):
    """
    modify atom positions by randomly choosing a H atom, 
    generating a random vector with length distributed according
    to normal distribution with centre mu and stdev sigma and
    perturbing that H atom's coordinates

    Parameters 
    ----------
    atoms: ase.Atoms obj 
           atoms object to modify coordinates of (needs to already have
           nuclear charges and coordinates set)
    mu: float 
        centre of normal distribution
    sigma: float 
           stdev of normal distribution

    Returns
    -------
    atoms: ase.Atoms obj
           atoms obj with coordinates of a random H atom modified
    """
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

def optimise_geometry(atoms, maxiter=1000, RT=0.1, savegif=None):
    """
    using PM7 in MOPAC

    Parameters
    ----------
    atoms: ase.Atoms object 
           Initial geometry to optimise 
    maxiter: int 
             number of iterations of MC optimisation
    RT: float 
        temperature for the Boltzmann distribution
    savegif: str 
             xyz file to save optimisation steps to

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

        # save to file
        if savegif:
            ase_tools.write_atoms_to_file(atoms_new, savegif)

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

    acceptance_ratio = (accept / maxiter) * 100
    print("acceptance ratio", acceptance_ratio, "over", maxiter, "iterations")

    # dump outfiles
    ase_tools.del_outfiles()

    return atoms_new, E_new


if __name__ == "__main__":
    xyz_in, n_steps, xyz_out, savegif = parse_args()
    atoms = ase_tools.init_atoms_from_file(xyz_in)
    atoms_new, E_new = optimise_geometry(atoms, maxiter=n_steps, savegif=savegif)
    print("final energy", E_new)
    ase_tools.write_atoms_to_file(atoms_new, xyz_out)

