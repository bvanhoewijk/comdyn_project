

# Cutoff threholds

## 5eqi -> 4zw9

### Setup files for MD experiment:


Aim is to generate a bunch of ITP files with thresholds. This script populates the 02_angstrom_benchmark folder and overwrites existing itp files if they are there.

```bash
PROJECT=/home/bas/projects/comdyn_project
cd ${PROJECT}
```


### Baseline
```bash
cp -r data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/
```
```bash
cp ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/5eqi_proa.itp
```


### 1 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/5eqi_proa.itp \
--nm 0.1 \
--martini1 --write
```

### 0.1 Angstrom threshold

```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/5eqi_proa.itp \
--nm 0.01 \
--martini1 --write
```

### 0.5 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/5eqi_proa.itp \
--nm 0.05 \
--martini1 --write
```

### 0.25 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/5eqi_proa.itp \
--nm 0.025 \
--martini1 --write
```

