#!/bin/bash

module load 2022
module load GROMACS/2021.6-foss-2022a

for i in 1 2 3 4 5
do
printf "1\n1\n" | gmx gyrate \
    -s rep${i}/step7_production_rep${i}.tpr \
    -f rep${i}/step7_production_rep${i}.xtc \
    -o gyration_rep${i}.xvg
done
