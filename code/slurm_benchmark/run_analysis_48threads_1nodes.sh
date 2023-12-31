#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --partition=genoa
#SBATCH --time=00:10:00
#SBATCH --exclusive
 
module load 2022
module load GROMACS/2021.6-foss-2022a

# Minimization
setenv GMX_MAXCONSTRWARN -1
export OMP_NUM_THREADS=1
# step6.0 - soft-core minimization
# If you encountered "There are 1 perturbed non-bonded pair interaction ......" error message,
# please modify rvdw and rcoulomb values from 1.1 to 2.0 in the step6.0_minimization.mdp file
#srun gmx grompp -f step7_production.mdp -o step7_production.tpr -c step6.6_equilibration.gro -p system.top -n index.ndx
#srun gmx_mpi mdrun -deffnm step7_production -pin on -g one_nodes_48threads.log 
mpirun gmx_mpi mdrun -deffnm step7_production -pin on -g one_nodes_48threads.log 
