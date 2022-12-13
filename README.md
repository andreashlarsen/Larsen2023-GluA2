# Larsen2023-GluA2
Scripts for paper: Larsen, Perozzo, Biggin, Bowie and Kastrup 2023 (in prep).

## Folder: PMF
PMF folder contains all scripts for free energy calculations usign potential of mean force (PMF).    
All-atomistic simulations.  
Workflow can be found in PMF/bash_scripts/RUN.sh, which is the master script that calls the other scripts.    
(Flow_UMBRELLA.sh is a script for running umbrella simulations on an HPC).    

## Folder: MetaD_BME
MetaD_BME folder contains all scripts for running coarse-grained simulations enhanced with metadynamics (MetaD).    
Trajectory was reweighted using the Bayesian/Maximum Entropy algorithm (BME) to get better consistency with small-angle neutron scattering (SANS) data.    
Workflow can be found in MetaD_BME/bash_scripts/Flow.sh, which is the master script that calls the other scripts.     

## Folder: Modeller    
Files for adding missing residues to PDB 5VHZ (dezensitized GluA2). Run using bash script script_modeller.sh.     

## Reproducibility and disclaimer
The scripts will not work "out of the box". Paths must be adjusted, file names and paths updated to the local system (and HPC), and programs installed (e.g. Python and GROMACS with PLUMED - see manuscript). But with some adjustment, the scripts should be sufficient to reproduce the results (disregarding differences due to inherent stochastic processes). 
