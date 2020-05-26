#!/usr/bin/env python3
"""
test case for the optimiser
"""
import tools
import optHs

def ethene_distorted():
    """
    this should be quite easy to optimise as 
    I only slightly distorted one of the C and 
    H positions in the molecule
    """
    file = "../data/dsgdb9nsd_000004/sample.xyz"
    OBMol = tools.read_OBMol(file)
    optHs.optimise_geometry(OBMol, maxiter=20)

if __name__ == "__main__":
    ethene_distorted()
