printf "0\n" | gmx trjconv -s step7_production.tpr -f step7_production.trr -o md_01_noPBC.xtc -pbc mol -ur compact

printf "4\n" | gmx trjconv -s step7_production.tpr -f md_01_noPBC.xtc -o md_01_noWAT.xtc -n index

printf "1 1\n" | gmx covar -s step7_production.tpr -f md_01_noWAT.xtc

printf "1 1\n" | gmx anaeig -s step7_production.tpr -f md_01_noWAT.xtc -2d -first 1 -last 2 -filt movie.pdb