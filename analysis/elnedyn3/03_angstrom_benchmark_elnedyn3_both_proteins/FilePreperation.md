# Unpacking zip files and folder prep

## 5eqi setup
```bash
cd 5eqi
tar -xzvf 5eqi_minimized_wildtype.tgz
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline
```

```bash
cd charmm-gui-0955694022/gromacs/
cp -r * ../../0.1nm/
cp -r * ../../0.05nm/
cp -r * ../../0.025nm/
cp -r * ../../0.01nm/
cp -r * ../../baseline/
```

## 4zw9 setup
```bash
cd 4zw9
tar -xzvf 4zw9_minimized_wildtype.tgz
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline
```

Copy files
```bash
cd charmm-gui-0990339221/gromacs/
cp -r * ../../0.1nm/
cp -r * ../../0.05nm/
cp -r * ../../0.025nm/
cp -r * ../../0.01nm/
cp -r * ../../baseline/
```

## 4zwc setup
```bash
cd 4zwc
tar -xzvf 4zwc_membrane.tgz
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline
```

Copy files
```bash
cd charmm-gui-0990049626/gromacs/
cp -r * ../../0.1nm/
cp -r * ../../0.05nm/
cp -r * ../../0.025nm/
cp -r * ../../0.01nm/
cp -r * ../../baseline/
```

## 4zwb setup
```bash
cd 4zwb
tar -xzvf 4zw9_minimized_wildtype.tgz
mkdir -p 0.01nm 0.025nm 0.05nm 0.1nm baseline
```

Copy files
```bash
cd charmm-gui-0989016529/gromacs/
cp -r * ../../0.1nm/
cp -r * ../../0.05nm/
cp -r * ../../0.025nm/
cp -r * ../../0.01nm/
cp -r * ../../baseline/
```

# ComDYN

cd /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both
```
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both
CODE=/home/bas/projects/comdyn_project/code

cd ${BASEPATH}
```

## 4ZW9
Do COMDYN for the 4ZW9 protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.025 0.05 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/baseline/4zw9_minimized_proa.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/baseline/5eqi_minimized_proa.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/${value}nm/4zw9_minimized_proa.itp \
    --nm ${value} \
    --write
done
```
## 4ZWC
Do COMDYN for the 4ZW9 protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.025 0.05 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwc/baseline/4zwc_proa.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/baseline/5eqi_minimized_proa.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwc/${value}nm/4zwc_proa.itp \
    --nm ${value} \
    --write
done
```

## 4ZWB
Do COMDYN for the 4ZWB protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.025 0.05 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwb/baseline/4zwb_proa.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/baseline/5eqi_minimized_proa.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwb/${value}nm/4zwb_proa.itp \
    --nm ${value} \
    --write
done
```


## 5EQI
Do COMDYN for the 5EQI protein for the four thresholds. Note: existing ITP files get backed-up.

```bash
for value in 0.1 0.025 0.05 0.01
do
    python ${CODE}/comdyn.py \
    --itp_file1 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/baseline/5eqi_minimized_proa.itp \
    --itp_file2 /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/baseline/step6_proa.itp \
    --itp_out /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/${value}nm/5eqi_minimized_proa.itp \
    --nm ${value} \
    --write
done
```

# Copy configuration

```bash
CODE=/home/bas/projects/comdyn_project/code
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    # MDP files
    cp ${CODE}/mdp/* /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/${value}/
    cp ${CODE}/mdp/* /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/${value}/
    cp ${CODE}/mdp/* /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwc/${value}/

    # Calibration script:
    cp ${CODE}/run_calibration.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/${value}/
    cp ${CODE}/run_calibration.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/${value}/
    cp ${CODE}/run_calibration.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwc/${value}/

    # Production run script:
    cp ${CODE}/slurm_start_simulation8h_propermpi_3reps.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zw9/${value}/
    cp ${CODE}/slurm_start_simulation8h_propermpi_3reps.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/4zwc/${value}/
    cp ${CODE}/slurm_start_simulation8h_propermpi_3reps.sh /home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both/5eqi/${value}/
done
```

# Run calibration

```bash
conda activate md_production
source /usr/local/gromacs/bin/GMXRC

BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both

# 4ZW9
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    cd ${BASEPATH}/4zw9/${value}
    bash run_calibration.sh
done

# 4ZWB
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    cd ${BASEPATH}/4zwb/${value}
    bash run_calibration.sh
done

# 5EQI:
for value in 0.1nm 0.025nm 0.05nm 0.01nm baseline
do
    cd ${BASEPATH}/5eqi/${value}
    bash run_calibration.sh
done
```

# Add stuff to git after calibration

## 4ZW9
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both

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
BASEPATH=/home/bas/projects/comdyn_project/analysis/02_angstrom_benchmark_both

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