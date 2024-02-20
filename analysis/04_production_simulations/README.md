# 1. Do COMDYN

The unfiltered itp files can be found in zip file. "all_membrane_systems_clean.zip". 

Note: The modified MDP files are also in the zip. These are not the default ones from CHARMMGUI. (`s/berendsen/[cv]rescale/`)

Unpack and run the following lines of code to remove the constaints. 

```bash
# Mutant 1
python ../../code/comdyn.py --itp_file1 4zw9_331trp/4zw9_trp331_proa.itp --itp_file2 5eqi_333trp/5eqi_trp333_proa.itp --write
mv output.itp 4zw9_331trp/4zw9_trp331_proa.itp 

python ../../code/comdyn.py --itp_file1 5eqi_333trp/5eqi_trp333_proa.itp --itp_file2 4zw9_331trp/4zw9_trp331_proa.itp --write
mv output.itp 5eqi_333trp/5eqi_trp333_proa.itp

# Mutant 2
python ../../code/comdyn.py --itp_file1 4zw9_89asp/4zw9_asp89_proa.itp --itp_file2 5eqi_91asp/5eqi_asp91_proa.itp --write
mv output.itp 4zw9_89asp/4zw9_asp89_proa.itp

python ../../code/comdyn.py --itp_file1 5eqi_91asp/5eqi_asp91_proa.itp --itp_file2 4zw9_89asp/4zw9_asp89_proa.itp --write
mv output.itp 5eqi_91asp/5eqi_asp91_proa.itp

######## Wildtype
# To keep: 2034/2545 79.9%
# To skip: 511/2545 20.1%
python ../../code/comdyn.py --itp_file1 4zw9_wildtype/4zw9_proa.itp --itp_file2 5eqi_wildtype/5eqi_proa.itp
mv output.itp 4zw9_wildtype/4zw9_proa.itp

# To keep: 2034/2318 87.7%
# To skip: 284/2318 12.3%
python ../../code/comdyn.py --itp_file1 5eqi_wildtype/5eqi_proa.itp --itp_file2 4zw9_wildtype/4zw9_proa.itp
mv output.itp 5eqi_wildtype/5eqi_proa.itp
```

# 2. Do equilibration

First copy the calibration script:
```bash
cp ../../code/run_calibration.sh 4zw9_331trp/
cp ../../code/run_calibration.sh 5eqi_333trp/
cp ../../code/run_calibration.sh 4zw9_89asp/
cp ../../code/run_calibration.sh 5eqi_91asp/
cp ../../code/run_calibration.sh 4zw9_wildtype/
cp ../../code/run_calibration.sh 5eqi_wildtype/
```

Now run:
```bash
current=`pwd`
source /usr/local/gromacs/bin/GMXRC
for protein in 4zw9_331trp  4zw9_89asp  4zw9_wildtype  5eqi_333trp  5eqi_91asp  5eqi_wildtype
do
cd ${protein}
bash ./run_calibration.sh
cd ${current}
done
```