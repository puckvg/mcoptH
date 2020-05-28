# Monte Carlo optimisation of H-atom positions

This library optimises the positions of H-atoms in a molecule, keeping the positions of the other atoms fixed. 
There is an option to use ASE's LBFGS optimiser at the end.

Intended usage: 
```
./optHs.py file_to_opt n_steps out_file save_gif extra_opt
```
where file_to_opt is the xyz file to optimise, n_steps is the number of steps taken by the optimiser, out_file is the 
name of the xyz outfile, save_gif is an option to save the gif of the MC optimisation, extra_opt is an option to add 
on an LBFGS optimisation routine at the end. For example:
```
./optHs.py ../data/dsgdb9nsd_000004/path_0_frag_0_parent.xyz 100 "opt.xyz" "n" "n"
```
takes 100 steps to optimise a file provided in data and saves it as opt.xyz.

```
./optHs.py ../data/dsgdb9nsd_000004/path_0_frag_0_parent.xyz 100 "opt.xyz" "n" "y"
```
takes 100 steps of the MC optimiser, then adds an LBFGS cycle and saves the final structure as opt.xyz.

## Installation and Dependencies 
This is dependent on: 
- numpy (do I need to tell you how to download numpy?)
- ase (see [their documentation](https://wiki.fysik.dtu.dk/ase/install.html) for installation instructions)



