# Monte Carlo optimisation of H-atom positions

This library optimises the positions of H-atoms in a molecule, keeping the positions of the other atoms fixed. 

Intended usage: 
```
./optHs.py file_to_opt n_steps out_file
```
where file_to_opt is the xyz file to optimise, n_steps is the number of steps taken by the optimiser, and 
out_file is the name of the xyz outfile. For example:
```
./optHs.py ../data/dsgdb9nsd_000004/path_0_frag_0_parent.xyz 100 "opt.xyz"
```
takes 100 steps to optimise a file provided in data and saves it as opt.xyz

## Installation and Dependencies 
This is dependent on: 
- numpy (do I need to tell you how to download numpy?)
- ase (see [their documentation](https://wiki.fysik.dtu.dk/ase/install.html) for installation instructions)

## TODO
- Seems good for poor initial guesses, but not good for good initial guesses 
- Add final geometry opt on top


