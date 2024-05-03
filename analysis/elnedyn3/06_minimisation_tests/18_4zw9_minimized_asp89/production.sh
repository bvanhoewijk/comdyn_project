#!/bin/bash
source /usr/local/gromacs/bin/GMXRC
gmx grompp -f step7_production.mdp -o step7_production.tpr -c step6.6_equilibration.gro -p system.top -n index.ndx
gmx mdrun -deffnm step7_production

