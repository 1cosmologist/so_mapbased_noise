#!/bin/bash

#SBATCH -A mp107b
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --time=00:30:00
#SBATCH --nodes=4
#SBATCH --job-name=lat_noise_gen
#SBATCH -o /pscratch/sd/s/shamikg/so_mapbased_noise/output/slurm_logs/%x_job%j.out

# Parallelization                                                                                                                                                                                        
export OMP_NUM_THREADS=16
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export NUMBA_NUM_THREADS=$OMP_NUM_THREADS

export JAX_PLATFORMS=cpu

let nnode=$SLURM_JOB_NUM_NODES
# 128 cores, 256 hardware threads                                                                                                                                                                        
let ntask_node=256/$OMP_NUM_THREADS
let ntask=$nnode*$ntask_node
let ncore=$OMP_NUM_THREADS

srun -N $nnode -n $ntask -c $ncore --cpu_bind=cores python /pscratch/sd/s/shamikg/so_mapbased_noise/code/run_generate_noise_lat.py

