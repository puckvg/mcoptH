#!/usr/bin/env python3
"""
Generate points on the surface of a sphere
"""
import numpy as np

def sample_sphere(r):
    """
    In spherical polar coordinates we have: 
        r: radius of sphere 
        theta: polar angle between [0, pi]
        phi: azimuthal angle between [0, 2pi]

    However if we uniformly sample we will get points
    densely distributed on the poles of the sphere. 
    To get points uniformly distributed on the sphere 
    surface, sample theta from cos(theta) where cos(theta)
    in [-1,1]
    """
    costheta = np.random.uniform(-1, 1)
    theta = np.arccos(costheta)
    phi = np.random.uniform(0, 2*np.pi)
    return np.array([r, theta, phi])

def convert_spherical_cartesian(spherical_coords):
    """
    Convert from spherical coordinates r, theta, phi
    where theta is the polar angle and phi is the 
    azimuthal angle to Cartesian coordinates x, y, z
    """
    r = spherical_coords[0]
    theta = spherical_coords[1]
    phi = spherical_coords[2]

    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    return np.array([x, y, z])

def sample_sphere_cartesian(r):
    sph_coords = sample_sphere(r)
    cartesian_coords = convert_spherical_cartesian(sph_coords)
    return cartesian_coords

if __name__ == "__main__":
    coords = sample_sphere_cartesian(0.9)
    print(coords)

