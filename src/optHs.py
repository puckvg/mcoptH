#!/usr/bin/env python3
"""
program to optimise positions of H atoms in a molecule
as read from an xyz file 
"""
import openbabel as ob
import numpy as np

import tools
import sphere_point_picking

def get_random(low=0, high=1, size=1):
    return np.random.uniform(low=low, high=high, size=size)

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
    # get random length 
    r = np.random.uniform()
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

def update_mol_coords(OBMol):
    atoms = ob.OBMolAtomIter(OBMol)
    nuclear_charges, coordinates = tools.get_nuclear_charges_coordinates(OBMol)
    H_indices = np.where(nuclear_charges == 1)[0]
    
    # change H coordinates
    coordinates[H_indices] = update_coords(coordinates[H_indices])

    return coordinates

def make_updated_mol(OBMol):
    """
    make a new OBMol with updated coordinates
    """
    new_mol = ob.OBMol()

    new_coordinates = update_mol_coords(OBMol)
    atoms = ob.OBMolAtomIter(OBMol)
    nuclear_charges = [atom.GetAtomicNum() for atom in atoms]

    for i in range(nuclear_charges):
        atom = mol.NewAtom()
        nuclear_charge = nuclear_charges[i]
        atom.SetAtomicNum(nuclear_charge)
        coords = new_coordinates[i]
        atom.SetVector(coords[0], coords[1], coords[2])

    assert len(ob.OBMolAtomIter(OBMol)) == len(ob.OBMolAtomIter(new_mol)), \
            "new mol is not the same size as the old one"
    return new_mol

def calc_energy(OBMol, method="UFF"):
    ff = ob.OBForceField.FindForceField(method)
    ff.Setup(OBMol)
    return ff.Energy()


def minimise_energy(OBMol, method="UFF", maxiter=1000, RT=0.2):
    # initialise variables
    mol_old = OBMol 
    E_old = 0.
    E_new = 0.
    dE = 0.
    # counter
    accept = 0

    # begin simulation 
    for i in range(maxiter):
        print("iteration ", i)
        # translation move 
        mol_new = make_updated_mol(mol_old)

        # calc energy of new geometry 
        E_new = calc_energy(mol_new)
        print("energy ", E_new)
        dE = E_new - E_old
        print("dE ", dE)

        # Metropolis-Hastings acceptance criterion
        boltzmann = get_boltzmann(dE, RT)
        if boltzmann >= get_random():
            # accept new geometry
            print("accept")
            mol_old = mol_new 
            E_old = E_new
            accept += 1

        # otherwise leave as is
        else:
            print("reject new geometry")

    acceptance_ratio = (accept / maxiter) * 100
    print("acceptance ratio ", acceptance_ratio)

    return mol_new 


if __name__ == "__main__":
    trial = "../data/dsgdb9nsd_000004/path_0_frag_0_parent.xyz"
    OBMol = tools.read_OBMol(trial)
    update_mol_coords(OBMol)
    #minimise_energy(OBMol)








