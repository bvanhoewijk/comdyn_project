#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=staging 
#SBATCH --time=01:00:00       # estimated runtime

module load 2023
module load GROMACS/2023.3-foss-2023a

for i in 1 2 3
do
printf "0" | gmx trjconv -s rep${i}/step7_production_rep${i}.tpr \
	-f rep${i}/step7_production_rep${i}.xtc \
	-o rep${i}/md_rep${i}_noPBC.xtc -pbc mol -ur compact

printf "0" | gmx trjconv -s rep${i}/step7_production_rep${i}.tpr \
	-f rep${i}/md_rep${i}_noPBC.xtc -n index -o rep${i}/md_rep${i}_noWAT.xtc \

# Full protein:
printf "1 1\n" | gmx covar -s rep${i}/step7_production_rep${i}.tpr \
	-f rep${i}/md_rep${i}_noWAT.xtc -o rep${i}/eigenval.xvg

# Only the TM part:
printf "21 21\n" | gmx covar -s rep${i}/step7_production_rep${i}.tpr \
	-f rep${i}/md_rep${i}_noWAT.xtc -o rep${i}/eigenval.tm.xvg -n index.bb


# Full TM helices:
printf "22 22\n" | gmx covar -s rep${i}/step7_production_rep${i}.tpr \
	        -f rep${i}/md_rep${i}_noWAT.xtc -o rep${i}/eigenval.alltm.xvg -n index.bb

# Full protein
printf "1 1\n" | gmx anaeig -s rep${i}/step7_production_rep${i}.tpr -f rep${i}/md_rep${i}_noWAT.xtc \
	-first 1 -last 2 -eig rep${i}/eigenval.xvg -2d 2dproj_rep${i}.xvg

# Only the TM AA:
printf "21 21\n" | gmx anaeig -s rep${i}/step7_production_rep${i}.tpr -f rep${i}/md_rep${i}_noWAT.xtc \
	-first 1 -last 2 -eig rep${i}/eigenval.tm.xvg -2d 2dproj_rep${i}.tm.xvg -n index.bb

# Full TM:
printf "22 22\n" | gmx anaeig -s rep${i}/step7_production_rep${i}.tpr -f rep${i}/md_rep${i}_noWAT.xtc \
	        -first 1 -last 2 -eig rep${i}/eigenval.alltm.xvg -2d 2dproj_rep${i}.alltm.xvg -n index.bb


done
