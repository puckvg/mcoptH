# Monte Carlo optimisation of H-atom positions

This library optimises the positions of H-atoms in a molecule, keeping the positions of the other atoms fixed. 
There is an option to use ASE's LBFGS optimiser at the end.

Intended usage: 
```
./optHs.py file_to_opt n_steps out_file save_gif extra_opt
```
where file_to_opt is the xyz file to optimise, n_steps is the number of steps taken by the optimiser, out_file is the 
name of the xyz outfile, save_gif is an option to save the xyz file of the MC optimisation to make a gif, extra_opt is an option to add 
on a steepest descent optimisation routine at the end. 

I suggest around 250 optimisation steps of the MC routine to get a decent initial guess for the steepest descent algorithm 
at the end.

## Installation and Dependencies 
This is dependent on: 
- numpy (do I need to tell you how to download numpy?)
- ase (see [their documentation](https://wiki.fysik.dtu.dk/ase/install.html) for installation instructions)
- a MOPAC installation to use within ase 

