echo 0 | gmx trjconv -s step7_production.tpr -f step7_production.xtc -o md_01_noPBC.xtc -pbc mol -ur compact
# Remove water:
echo 0 | gmx trjconv -s step7_production.tpr -f md_01_noPBC.xtc -o md_01_noWAT.xtc -n index
# Covariance
printf "1 1\n" | gmx covar -s step7_production.tpr -f md_01_noWAT.xtc
# Map points on eigenvectors/values:
printf "1 1\n" | gmx anaeig -s step7_production.tpr -f md_01_noWAT.xtc -2d -first 1 -last 2 -filt movie.pdb -skip 1
