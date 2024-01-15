# 1 Folder setup

Folders:
1. MD: Results of molecular dynamics runs
2. MD_analysis: Notebooks + generated files for the purpose of interpretation of the data. 

Set correct directory:
```bash
PROJECT=/home/bas/projects/comdyn_project
# PROJECT=/home/bvanhoewijk/project/comdyn_project
cd ${PROJECT}
```


Copy the data:
```bash
cp -r ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/
cp -r ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/
cp -r ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/
cp -r ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/
cp -r ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/
```

# 2 Setup 

## 1 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/5eqi_proa.itp \
--nm 0.1 \
--martini1 --write
```

## 0.5 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/5eqi_proa.itp \
--nm 0.05 \
--martini1 --write
```

## 0.25 Angstrom threshold
```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/5eqi_proa.itp \
--nm 0.025 \
--martini1 --write
```

## 0.1 Angstrom threshold

```bash
python code/comdyn.py \
--itp_file1 ${PROJECT}/data/charmmgui_pdb/5eqi_membrane/gromacs/5eqi_proa.itp \
--itp_file2 ${PROJECT}/data/charmmgui_pdb/4zw9_membrane/gromacs/4zw9_proa.itp \
--itp_out ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/5eqi_proa.itp \
--nm 0.01 \
--martini1 --write
```


# 3 Run Energy minimization and equilibration

## 3.1 Copy required files
Copy calibration script:
```bash
cp ${PROJECT}/code/run_calibration.sh ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/
cp ${PROJECT}/code/run_calibration.sh ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/
cp ${PROJECT}/code/run_calibration.sh ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/
cp ${PROJECT}/code/run_calibration.sh ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/
cp ${PROJECT}/code/run_calibration.sh ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/
```

Overwrite MDP files:
```bash
cp ${PROJECT}/code/mdp/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/
cp ${PROJECT}/code/mdp/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/
cp ${PROJECT}/code/mdp/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/
cp ${PROJECT}/code/mdp/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/
cp ${PROJECT}/code/mdp/* ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/
```

## 3.2 Run

### Baseline
```bash
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline/
bash run_calibration.sh
```


### 1.0 angstrom
```bash
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/
bash run_calibration.sh
```

### 0.5 angstrom
```bash
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/
bash run_calibration.sh
```

### 0.25 angstrom
```bash
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/
bash run_calibration.sh
```

### Everything:
```bash
PROJECT=/home/bas/projects/comdyn_project
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/ && bash run_calibration.sh
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/ && bash run_calibration.sh
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/ && bash run_calibration.sh
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/ && bash run_calibration.sh
```


# 4 Add relevant items to git and push

```bash
cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.1_angstrom/
git add step6.6_equilibration.gro index.ndx step7_production.mdp system.top 5eqi_proa.itp toppar

cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.25_angstrom/
git add step6.6_equilibration.gro index.ndx step7_production.mdp system.top 5eqi_proa.itp toppar

cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/0.5_angstrom/
git add step6.6_equilibration.gro index.ndx step7_production.mdp system.top 5eqi_proa.itp toppar

cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/1.0_angstrom/
git add step6.6_equilibration.gro index.ndx step7_production.mdp system.top 5eqi_proa.itp toppar

cd ${PROJECT}/analysis/02_angstrom_benchmark/MD/baseline
git add step6.6_equilibration.gro index.ndx step7_production.mdp system.top 5eqi_proa.itp toppar

git commit -m 'calibration'
```
