#!/bin/bash

# Minimization
export GMX_MAXCONSTRWARN=-1
# step6.0 - soft-core minimization
# If you encountered "There are 1 perturbed non-bonded pair interaction ......" error message, 
# please modify rvdw and rcoulomb values from 1.1 to 2.0 in the step6.0_minimization.mdp file
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c step5_charmm2gmx.pdb -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
gmx mdrun -deffnm step6.0_minimization

# # step6.1
gmx grompp -f step6.1_minimization.mdp -o step6.1_minimization.tpr -c step6.0_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx -maxwarn 1
gmx mdrun -deffnm step6.1_minimization
unset GMX_MAXCONSTRWARN

# Equilibration
cnt=2
cntmax=6

while [ $cnt -le $cntmax ]; do
    pcnt=$((cnt - 1))
    if [ $cnt -eq 2 ]; then
        gmx grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_minimization.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx
    else
        gmx grompp -f step6.${cnt}_equilibration.mdp -o step6.${cnt}_equilibration.tpr -c step6.${pcnt}_equilibration.gro -r step5_charmm2gmx.pdb -p system.top -n index.ndx
    fi
    gmx mdrun -deffnm step6.${cnt}_equilibration
    ((cnt++))
done
