# Larsen2022-GluA2
Scripts for paper: Larsen, Biggin, Jensen 2022

## Folder: PMF
Folder contains all scripts for free energy calculations usign potential of mean force (PMF).    
All-atomistic simulations.  
Workflow in RUN.sh, which is the master script that calls the other scripts.   
Flow_UMBRELLA.sh is a script for running umbrella simulations on a HPC.   

## Folder: MetaD_BMEMetaD_BME
Folder contains all scripts for running coarse-grained simulations enhanced with metadynamics (MetaD).    
Trajectory was reweighted using the Bayesian/Maximum Entropy algorithm (BME) to get better consistency with small-angle neutron scattering (SANS) data.    

Workflow in MetaD_BME/bash_scripts/Flow.sh, which is the master script that calls the other scripts.     

## Reproducibility and disclaimer
scripts will not work "out of the box". Paths must be adjusted, file names and paths updated to the local system, and programs installed (e.g. Python and GROMACS with PLUMED - see manuscript). But with some adjustment, the scripts should be sufficient to reproduce the results (disregarding differences due to inherent stochastic processes). 
