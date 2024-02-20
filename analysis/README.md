
```bash
PROJECT=/home/bvanhoewijk/project/comdyn_project/analysis/02_angstrom_benchmark/MD/
for item in 0.1_angstrom 0.25_angstrom 0.5_angstrom 1.0_angstrom baseline
do
cd ${PROJECT}${item}
bash rmsd.sh
bash distance.sh
done
```
