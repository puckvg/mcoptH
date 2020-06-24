# Markov-Chain Monte Carlo (MCMC) optimisation of H-atom positions

## Theory
MCMC methods are a class of algorithms for sampling from a probability distribution. 
Markov chains suggest that the future state of aystem depends on the current state only, rather than the entire history. 

### Metropolis-Hastings algorithm
Assume the Markov chain is in some state $X_n = i$. We generate $X_{n+1}$ using the following algorithm: 
1. Choose a proposed state $j$ 
2. Compute the acceptance probability $\alpha_{ij} = \mathrm{min}(1, \exp^{-\Beta (E_j - E_i)})$ where $\Beta = \frac{1}{RT}$ for molecules
3. Generate a uniform random number $U \sim \mathrm{Uniform(0,1)}$. If $U < \alpha_{ij}$, accept the move and set $X_{n+1} = j$. Otherwise reject the move and keep $X_{n+1}=X_n$

## This library
This library optimises the positions of H-atoms in a molecule using the Metropolis-Hastings algorithm described above, keeping the positions of the other atoms fixed. There is an option to use ASE's LBFGS optimiser at the end.

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

## TODO
Sometimes get this error "ase.calculators.calculator.PropertyNotImplementedError: forces not present in this calculation
"

