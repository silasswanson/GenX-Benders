#!/bin/bash

#SBATCH --job-name=TDR_sweep_ce
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ss9746@princeton.edu
#SBATCH --partition=cpu

module purge
module load julia/1.11.1
module load gurobi/10.0.1

cd /scratch/gpfs/JENKINS/ss9746/GenX-Benders

julia --project=. CE_TDR_RP_sweep.jl
