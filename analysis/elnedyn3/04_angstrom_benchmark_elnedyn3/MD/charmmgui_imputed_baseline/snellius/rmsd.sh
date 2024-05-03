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
done
