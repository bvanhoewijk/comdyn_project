#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=genoa
#SBATCH --time=1:00:00

module load 2022
module load GROMACS/2021.6-foss-2022a

for rep in 1 2 3 4 5
do
srun printf "1\n" | gmx trjconv \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f rep${rep}/step7_production_rep${rep}.xtc \
	-o md_01_noPBC.rep${rep}.xtc -pbc mol -ur compact

srun printf "1\n" | gmx trjconv \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f md_01_noPBC.rep${rep}.xtc \
	-o md_01_noWAT.rep${rep}.xtc

srun printf "17\n17\n" | gmx covar \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f md_01_noWAT.rep${rep}.xtc \
	-o rep${rep}.eigenval.xvg -v rep${rep}.eigenvec.trr \
	-av rep${rep}.average.pdb -n index

srun printf "17\n17\n" | gmx anaeig \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f md_01_noWAT.rep${rep}.xtc \
	-v rep${rep}.eigenvec.trr \
	-eig rep${rep}.eigenval.xvg \
	-2d 2dproj_rep${rep}.xvg -first 1 -last 2 -n index
done
