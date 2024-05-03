#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=genoa
#SBATCH --time=12:00:00
module load 2023
module load GROMACS/2023.3-foss-2023a

#srun gmx grompp -f step7_production.mdp -o step7_production.tpr -c step6.6_equilibration.gro -p topol.top -n index.ndx
srun gmx mdrun -v -deffnm step7_production -nt 64 -pin on
