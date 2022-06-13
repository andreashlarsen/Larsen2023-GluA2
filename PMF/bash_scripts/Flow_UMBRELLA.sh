#!/bin/bash 

## Define if you are in testing mode (on short queue) or running

## Set the job name.
#SBATCH --job-name=GLUA2_FOLDER

## Set the number of nodes and cores. This shouldn't need changing.
## The maximum number of nodes which can be used on axon is 1.
#SBATCH --nodes=1

## Set the number of tasks on each node, this is usually the number of program executions done.
## For example, if you are running an mpi run, this would reflect the number passed to the -np flag.
#SBATCH --ntasks-per-node=1

## Set the number of cores per tasks. This usually reflects the number of threads (often openmp) that
## are being assigned per task. Benchmarking has shown that setting tasks to 1 and cpus-per-task to
## 16 for atomistic simulations and 8 for CG is a good starting point for getting maximum efficiency.
#SBATCH --cpus-per-task=6

## Set the number of GPUs to be used. In most cases this will be set to 1.
#SBATCH --gres=gpu:1

## IMPORTANT: set GPU binding, otherwise your jobs will clash with other users'
#SBATCH --gres-flags=enforce-binding

## Select the queues you will be running on (sansom: gpu-sansom,gpu-sansom2 biggin: gpu-biggin,gpu-biggin2) 
##SBATCH -p gpu-sm-short
#SBATCH -p gpu-sansom3,gpu-sansom4

## Select the max amount of time this job will run (48h for gpu-sansom, 3h for shor)
##SBATCH --time=3:00:00
#SBATCH --time=48:00:00

## Set an array of jobs you want to run
#SBATCH --array=1-NFRAMES%6

source /etc/profile.d/modules.sh

module purge
module load apps/gromacs/2018.6-plumed_2.4.4-GPU-KEPLER

## Note if running an energy minimisation, CG or using energy groups, you need to add the -ntmpi 1 for gromacs 2019 and above
#gmx mdrun -deffnm umbrella_step${SLURM_ARRAY_TASK_ID} -v -ntomp ${SLURM_CPUS_PER_TASK} -ntmpi 1 -quiet -nsteps 1000
gmx mdrun -deffnm umbrella_step${SLURM_ARRAY_TASK_ID} -v -ntomp ${SLURM_CPUS_PER_TASK} -ntmpi 1 -quiet
