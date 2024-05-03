# Unpacking zip files and folder prep

```bash
mkdir 5eqi
cd 5eqi
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline no_rubber all_common
cd ..

mkdir 4zw9
cd 4zw9
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline no_rubber all_common 
cd ..

```

## Unpack and copy files:
```bash
tar -xzvf 4zw9_elnedyn22.tgz
for item in 0.01nm  0.025nm  0.05nm  0.1nm  baseline
do
cp -r charmm-gui-1023580519/gromacs/* 4zw9/${item}/
done

rm -rf charmm-gui-1023580519
```

```bash
tar -xzvf 5eqi_elnedyn22.tgz
for item in 0.01nm  0.025nm  0.05nm  0.1nm  baseline
do
cp -r charmm-gui-1023500892/gromacs/* 5eqi/${item}/
done

rm -rf charmm-gui-1023500892
```

# ComDYN

```
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22
CODE=/home/bas/projects/comdyn_project/code

cd ${BASEPATH}
```

## 4ZW9
Do COMDYN for the 4ZW9 protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.05 0.025 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/${value}nm/PROA_P.itp \
    --nm ${value} \
    --write
done

for value in 0.1 0.05 0.025 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/${value}nm/PROA_P.itp \
    --nm ${value} >> 4zw9_log.txt
done 
```

## 5EQI
Do COMDYN for the 5EQI protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.05 0.025 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/${value}nm/PROA_P.itp \
    --nm ${value} \
    --write
done


for value in 0.1 0.05 0.025 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/${value}nm/PROA_P.itp \
    --nm ${value} >> 5eqi_log.txt
done
```

# ALL common
```bash
python ${CODE}/comdyn.py \
--itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
--itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
--itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/all_common/PROA_P.itp \
--nm 100 \
--write

python ${CODE}/comdyn.py \
--itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
--itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
--itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/all_common/PROA_P.itp \
--nm 100 \
--write
```

# Copy configuration

```bash
CODE=/home/bas/projects/comdyn_project/code
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline all_common
do
    # Config files
    cp ${BASEPATH}/config/* /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/${value}/
    cp ${BASEPATH}/config/* /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/${value}/
done
```

# Run calibration

```bash
conda activate md_production
source /usr/local/gromacs/bin/GMXRC

BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22

# 4ZW9
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    cd ${BASEPATH}/4zw9/${value}
    bash run_calibration.sh
done

# 5EQI:
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    cd ${BASEPATH}/5eqi/${value}
    bash run_calibration.sh
done
```


# All common
```bash
value=all_common
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22
for protein in 5eqi 4zw9
do
    cd ${BASEPATH}/${protein}/all_common
    bash run_calibration.sh
done
```



# Add calibrated datasets to git version control

## 4ZW9
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22

for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    git add ${BASEPATH}/4zw9/${value}/*.ndx
    git add ${BASEPATH}/4zw9/${value}/*.sh
    git add ${BASEPATH}/4zw9/${value}/*.itp
    git add ${BASEPATH}/4zw9/${value}/system.top
    git add ${BASEPATH}/4zw9/${value}/toppar/*
    git add -f ${BASEPATH}/4zw9/${value}/*.pdb
    git add -f ${BASEPATH}/4zw9/${value}/step6.6*.gro
    git add -f ${BASEPATH}/4zw9/${value}/step7_production.mdp
done
```

## 5EQI
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22

for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    git add ${BASEPATH}/5eqi/${value}/*.ndx
    git add ${BASEPATH}/5eqi/${value}/*.sh
    git add ${BASEPATH}/5eqi/${value}/*.itp
    git add ${BASEPATH}/5eqi/${value}/system.top
    git add ${BASEPATH}/5eqi/${value}/toppar/*
    git add -f ${BASEPATH}/5eqi/${value}/*.pdb
    git add -f ${BASEPATH}/5eqi/${value}/step6.6*.gro
    git add -f ${BASEPATH}/5eqi/${value}/step7_production.mdp
done
```


# Run benchmark on Snellius
Change to the 4zw9 directory:

```bash
cd 4zw9

for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
cd ${value} 
sbatch slurm_start_simulation8h_propermpi_3reps.sh
cd ..
done
```

Change to the 5eqi directory:

```bash
cd 5eqi

for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
cd ${value} 
sbatch slurm_start_simulation8h_propermpi_3reps.sh
cd ..
done
```

# Extra

```bash
CODE=/home/bas/projects/comdyn_project/code

for value in 0.15 0.2 
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/${value}nm/PROA_P.itp \
    --nm ${value} 
done

for value in 0.15 0.2 1.0
do
    python ${CODE}/comdyn.py \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/5eqi/baseline/PROA_P.itp \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_elnedyn22/4zw9/baseline/PROA_P.itp \
    --nm ${value} 
done | grep "To keep"
```