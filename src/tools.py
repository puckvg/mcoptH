#!/usr/bin/env python3 
"""
general tools I keep re-programming
"""

import numpy as np
import matplotlib.pyplot as plt 
import itertools

import os
import pickle 

from copy import deepcopy

import openbabel as ob

import qml
from qml.kernels import laplacian_kernel
from qml.kernels import kpca

from sklearn import decomposition

from nuclear_charge import ELEMENT


# file IO functions 
def strip_string(string):
    """
    strip string after last /
    """
    return string.split("/")[-1]

def strip_extension(file):
    """
    Get filename without extension
    """
    name = file.split(".")[0]
    return name 

def get_raw_filename(file):
    """
    Strip extension and beginning of path
    """
    filename = strip_string(file)
    return strip_extension(filename)

def get_file_lists(directory, sort=True):
    if sort:
        files = sorted(os.listdir(directory))
    else:
        files = os.listdir(directory)

    return [directory + file for file in files]

def get_N_files(directory, N, sort=True):
    files = get_file_lists(directory, sort=True)
    return files[:N]

def write_list_to_file(list_to_write, filename):
    with open(filename, "w") as f:
        for item in list_to_write:
            f.write(item)
            f.write("\n")
    return

def read_list_from_file(filename): 
    with open(filename, "r") as f:
        lines = f.readlines()

    lines = [line.strip() for line in lines]
    return lines

def pickle_load(file):
    return pickle.load(open(file, "rb"))

def pickle_dump(data, filename):
    if not filename.endswith(".pickle"):
        filename = filename + ".pickle"
    return pickle.dump(data, open(filename, "wb"))

def check_make_directory(path, exist_ok=True):
    """
    Check a directory exists. If it doesn't, make it.
    """
    os.makedirs(path, exist_ok)
    os.chmod(path, 0o777)
    return

def check_path_format(path):
    """
    Check path ends with a "/", if not, add one
    """
    if not path.endswith("/"):
        path = path + "/"
    return path 

# lists 
def get_index_from_list(list_to_search, entry):
    return list_to_search.index(entry)

def get_index_from_file(file_to_search, entry):
    lines = read_list_from_file(file_to_search)
    return get_index_from_list(lines, entry) 

def remove_item_from_list(list_to_modify, index):
    modified_list = deepcopy(list_to_modify)

    if type(index) == int:
        modified_list.pop(index)
    else:
        modified_list.remove(index)

    return modified_list

def init_list_of_lists(length):
    return [[] for i in range(length)]

# dictionary
def init_dict_from_list(list_keys):
    return {k:[] for k in list_keys}


# array functions
def return_min_size(*args):
    sizes = []
    for item in args: 
        sizes.append(len(item))

    min_size = min(sizes)
    reshaped = []
    for item in args: 
        item = item[:min_size]
        reshaped.append(item)

    return reshaped

def shuffle_list_indices(list_to_shuffle):
    random_indices = np.arange(len(list_to_shuffle))
    np.random.shuffle(random_indices)
    return random_indices

def concatenate_many_arrays(*args):
    a = np.concatenate((args[0], args[1]))
    for b in a[2:]:
        a = np.concatenate((a, b))
    return a

def flatten_first_dim(array):
    shape = array.shape
    first_dim = (shape[0]*shape[1],)
    other_dims = shape[2:]
    return np.reshape(array, first_dim + other_dims)

def l2_norm(a, b):
    return np.linalg.norm(a-b)
    
def return_matching_indices(arr, cond):
    return [idx[0] for idx in np.argwhere(arr == cond)]

def check_item_not_in_list(item, item_list, rtol=1e-3):
    """
    check item isn't already in list 
    """

    if np.any([np.allclose(item, other, rtol=rtol) for other in item_list]):
        return False 

    return True 

# mol functions
def get_heavy_atoms_coordinates(nuclear_charges, coordinates):
    heavy_indices = np.where(np.array(nuclear_charges) != 1)
    nuclear_charges = nuclear_charges[heavy_indices].ravel()
    coordinates = coordinates[heavy_indices]
    return nuclear_charges, coordinates

def write_xyz(nuclear_charges, coordinates, filename):
    n_atoms = len(nuclear_charges)
    if n_atoms == 0:
        return 
    atom_types = [ELEMENT[n] for n in nuclear_charges]

    with open(filename, "w") as f:
        f.write(str(n_atoms) + "\n\n")

        for atom, coords in zip(atom_types, coordinates):
            line = atom + "\t" + str(coords[0]) + " " +\
                   str(coords[1]) + " " + str(coords[2]) +\
                   "\n"
            f.write(line)

    print("file saved to ", filename)
    return

# ML functions
def split_data(data, n_train, n_test):
    train = data[:n_train]
    test = data[n_train:n_train+n_test]
    return train, test

def train_dump_model(filename, X, y, estimator):
    estimator.fit(X, y)
    pickle_dump(estimator, filename)
    return

def get_score_statistics(scores, axis=1):
    """
    scores have dimensions n_ticks, n_cv_folds 
    """
    mean = np.mean(scores, axis=axis)
    std = np.std(scores, axis=axis)
    return mean, std

def get_laplacian_kPCA(X, n=2, sigma=2e5):
    """
    kPCA with laplacian kernel
    """
    # TODO get max value of PCA 
    K = laplacian_kernel(X, X, sigma)
    pcas_qml = kpca(K, n)
    return pcas_qml

def get_PCA(X, n):
    pca = decomposition.PCA(n)
    pca.fit(X)
    X = pca.transform(X)
    return X

# OB functions 
def read_OBMol(filename, extension="xyz"):
    mol = ob.OBMol()

    obConversion = ob.OBConversion()
    obConversion.SetInFormat(extension)
    read_ok = obConversion.ReadFile(mol, filename)
    if not read_ok:
        raise Exception("Could not read file {}".format(filename))
        # or return None 

    return mol

def write_OBMol(OBMol, outfile, outformat="xyz"):
    obConversion = ob.OBConversion()
    obConversion.SetOutFormat(outformat)
    obConversion.WriteFile(OBMol, outfile)

    return

def get_nuclear_charges_coordinates(OBMol):
    atoms = ob.OBMolAtomIter(OBMol)
    nuclear_charges = []
    coordinates = []

    for atom in atoms:
        nuclear_charge = atom.GetAtomicNum()
        nuclear_charges.append(nuclear_charge)

        coord = get_coordinates_atom(atom)
        coordinates.append(coord)
    
    return np.array(nuclear_charges), np.array(coordinates)

def get_coordinates_atom(OBAtom):
    x = OBAtom.GetX()
    y = OBAtom.GetY()
    z = OBAtom.GetZ()
    coords = np.array([x, y, z])
    return coords 

def opt_OBMol(OBMol, forcefield="MMFF94", nsteps=1000):
    """
    Parameters 
    ----------
    OBMol: OB Mol 
           molecule to optimise
    forcefield: str
                type of forcefield 
    nsteps: int 
            number of steps to use in optimisation

    Returns 
    -------
    OBMol: OB Mol
           molecule with optimised geometry
    """
    FF = ob.OBForceField.FindType(forcefield)
    FF.Setup(OBMol)
    FF.ConjugateGradients(nsteps)
    FF.UpdateCoordinates(OBMol)

    return OBMol

def opt_geometry(nuclear_charges, coordinates, forcefield="UFF",
                 nsteps=1000):
    """
    hack to write xyz file with nuclear charges + coordinates 
    in order to geometry opt within python
    """
    # write tmp file
    write_xyz(nuclear_charges, coordinates, "tmp.xyz")
    OBMol = read_OBMol("tmp.xyz")

    # delete tmp file
    os.remove("tmp.xyz")
    print("deleted tmp file")

    # opt
    OBMol = opt_OBMol(OBMol)

    # get nuclear charges and coordinates 
    nuclear_charges, coordinates = get_nuclear_charges_coordinates(OBMol)
    return nuclear_charges, coordinates




