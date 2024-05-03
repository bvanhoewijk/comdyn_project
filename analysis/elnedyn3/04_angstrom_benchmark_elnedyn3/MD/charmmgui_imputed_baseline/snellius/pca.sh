#!/bin/bash

module load 2022
module load GROMACS/2021.6-foss-2022a

for i in 1 2 3 4 5
do
printf "17\n17\n" | gmx rms \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -o rmsd_rep${i}.xvg \
    -n index

printf "16\n16\n" | gmx covar \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f md_01_noWAT.rep${rep}.xtc \
	-o rep${rep}.eigenval.xvg -v rep${rep}.eigenvec.trr \
	-av rep${rep}.average.pdb -n index 2> /dev/null

printf "16\n16\n" | gmx anaeig \
	-s rep${rep}/step7_production_rep${rep}.tpr \
	-f md_01_noWAT.rep${rep}.xtc \
	-v rep${rep}.eigenvec.trr \
	-eig rep${rep}.eigenval.xvg \
	-2d 2dproj_rep${rep}.xvg -first 1 -last 2 -n index 2> /dev/null
done
