
# 1 Setup

Make 6 folders to store the files:
```
mkdir -p 4zw9 4zw9_asp89 4zw9_trp331
mkdir -p 5eqi 5eqi_asp91 5eqi_trp333
```

## 1.1 4ZW9 unpack and copy 
```bash
tar -xzvf 4ZW9_ASP89_elnedyn22.tgz
cp -r charmm-gui-1135530416/gromacs/* ./4zw9_asp89/
rm -rf charmm-gui-1135530416
```


```bash
tar -xzvf 4ZW9_TRP331_elnedyn22.tgz
cp -r charmm-gui-1135530095/gromacs/* ./4zw9_trp331/
rm -rf charmm-gui-1135530095
```

```bash
tar -xzvf 4zw9_elnedyn22.tgz 
cp -r charmm-gui-1023580519/gromacs/* ./4zw9/
rm -rf charmm-gui-1023580519
```

## 1.2 5EQI unpack and copy
```bash
tar -xzvf 5EQI_ASP91_elnedyn22.tgz
cp -r charmm-gui-1135480010/gromacs/* ./5eqi_asp91/
rm -rf charmm-gui-1135480010
```

```bash
tar -xzvf 5EQI_TRP333_elnedyn22.tgz
cp -r charmm-gui-1135479732/gromacs/* ./5eqi_trp333/
rm -rf charmm-gui-1135479732
```

```bash
tar -xzvf 5eqi_elnedyn22.tgz
cp -r charmm-gui-1023500892/gromacs/* ./5eqi/
rm -rf charmm-gui-1023500892
```

# 2 ComDYN

## 2.1 Wildtype
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22
CODE=/home/bas/projects/comdyn_project/code

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/4zw9/PROA_P.itp \
--itp_file2 ${BASEPATH}/5eqi/PROA_P.itp \
--itp_out ${BASEPATH}/4zw9/PROA_P.new.itp \
--nm 0.1 --write

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/5eqi/PROA_P.itp \
--itp_file2 ${BASEPATH}/4zw9/PROA_P.itp \
--itp_out ${BASEPATH}/5eqi/PROA_P.new.itp \
--nm 0.1 --write
```

## 2.2 Mutant1
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22
CODE=/home/bas/projects/comdyn_project/code

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/4zw9_asp89/PROA_P.itp \
--itp_file2 ${BASEPATH}/5eqi_asp91/PROA_P.itp \
--itp_out ${BASEPATH}/4zw9_asp89/PROA_P.new.itp \
--nm 0.1 --write

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/5eqi_asp91/PROA_P.itp \
--itp_file2 ${BASEPATH}/4zw9_asp89/PROA_P.itp \
--itp_out ${BASEPATH}/5eqi_asp91/PROA_P.new.itp \
--nm 0.1 --write
```

## 2.3 Mutant2
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22
CODE=/home/bas/projects/comdyn_project/code

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/4zw9_trp331/PROA_P.itp \
--itp_file2 ${BASEPATH}/5eqi_trp333/PROA_P.itp \
--itp_out ${BASEPATH}/4zw9_trp331/PROA_P.new.itp \
--nm 0.1 --write

python ${CODE}/comdyn.py \
--itp_file1 ${BASEPATH}/5eqi_trp333/PROA_P.itp \
--itp_file2 ${BASEPATH}/4zw9_trp331/PROA_P.itp \
--itp_out ${BASEPATH}/5eqi_trp333/PROA_P.new.itp \
--nm 0.1 --write
```

## 2.4
Update top files:
```bash
# sed -i is in place
find . -name "system.top" | xargs sed -i "s/PROA_P.itp/PROA_P.new.itp/"
```


# 3 Copy configuration files
```bash
CODE=/home/bas/projects/comdyn_project/code
BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22
for value in 4zw9 4zw9_asp89 4zw9_trp331 5eqi 5eqi_asp91 5eqi_trp333
do
    # Config files
    cp ${BASEPATH}/config/* ${BASEPATH}/${value}/
done
```

# 4 Calibration
```bash
conda activate md_production
source /usr/local/gromacs/bin/GMXRC

BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22

# DO THE THING:
for value in 4zw9 4zw9_asp89 4zw9_trp331 5eqi 5eqi_asp91 5eqi_trp333
do
    cd ${BASEPATH}/${value}/
    bash run_calibration.sh
done
```

# Save data to git
```bash
BASEPATH=/home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22

for value in 4zw9 4zw9_asp89 4zw9_trp331 5eqi 5eqi_asp91 5eqi_trp333
do
    git add ${BASEPATH}/${value}/*.ndx
    git add ${BASEPATH}/${value}/*.sh
    git add ${BASEPATH}/${value}/*.itp
    git add ${BASEPATH}/${value}/system.top
    git add ${BASEPATH}/${value}/toppar/*
    git add -f ${BASEPATH}/${value}/*.pdb
    git add -f ${BASEPATH}/${value}/step6.6*.gro
    git add -f ${BASEPATH}/${value}/step7_production.mdp
done
```