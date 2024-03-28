#!/bin/bash

module load 2023
module load GROMACS/2023.3-foss-2023a

for i in 1 2 3
do
printf "17\n" | gmx distance \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -oall distance1_rep${i}.xvg \
    -n index.bb

printf "18\n" | gmx distance \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -oall distance2_rep${i}.xvg \
    -n index.bb

printf "19\n" | gmx distance \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -oall distance3_rep${i}.xvg \
    -n index.bb

printf "20\n" | gmx distance \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -oall distance4_rep${i}.xvg \
    -n index.bb
done

