# 1. Do COMDYN

The unfiltered itp files can be found in zip file. "all_membrane_systems_clean.zip". 

Note: The modified MDP files are also in the zip. These are not the default ones from CHARMMGUI. (`s/berendsen/[cv]rescale/`)

Unpack and run the following lines of code to remove the constaints. 

## Mutant 1
To keep: 1908/2461 77.5%
To skip: 553/2461 22.5%

RB where AA not other file : 35
RB in both files. Size diff < threshold : 1908
RB in both files. Size diff >= threshold : 91
RB not in other file : 427

```bash
python ../../code/comdyn.py \
--itp_file1 ./4zw9_331trp/4zw9_minimized_trp331_proa.itp \
--itp_file2 ./5eqi_333trp/5eqi_minimized_trp333_proa.itp \
--itp_out ./4zw9_331trp/4zw9_minimized_trp331_proa.itp.new \
--write

python ../../code/comdyn.py \
--itp_file1 ./5eqi_333trp/5eqi_minimized_trp333_proa.itp \
--itp_file2 ./4zw9_331trp/4zw9_minimized_trp331_proa.itp \
--itp_out ./5eqi_333trp/5eqi_minimized_trp333_proa.itp.new \
--write


```

## Mutant 2
```bash
python ../../code/comdyn.py \
    --itp_file1 4zw9_89asp/4zw9_minimized_asp89_proa.itp \
    --itp_file2 5eqi_91asp/5eqi_minimized_asp91_proa.itp \
    --itp_out 4zw9_89asp/4zw9_minimized_asp89_proa.itp.new \
    --write

python ../../code/comdyn.py \
    --itp_file1 5eqi_91asp/5eqi_minimized_asp91_proa.itp \
    --itp_file2 4zw9_89asp/4zw9_minimized_asp89_proa.itp \
    --itp_out 5eqi_91asp/5eqi_minimized_asp91_proa.itp.new \
    --write


```


## Wildtype
```bash
python ../../code/comdyn.py \
    --itp_file1 4zw9_wildtype/step6_proa.itp \
    --itp_file2 5eqi_wildtype/5eqi_minimized_proa.itp \
    --itp_out 4zw9_wildtype/step6_proa.itp.new \
    --write

python ../../code/comdyn.py \
    --itp_file1 5eqi_wildtype/5eqi_minimized_proa.itp \
    --itp_file2 4zw9_wildtype/step6_proa.itp \
    --itp_out 5eqi_wildtype/5eqi_minimized_proa.itp.new \
    --write

mv ./4zw9_wildtype/step6_proa.itp.new 4zw9_wildtype/step6_proa.itp
mv ./5eqi_wildtype/5eqi_minimized_proa.itp.new 5eqi_wildtype/5eqi_minimized_proa.itp
mv ./4zw9_331trp/4zw9_minimized_trp331_proa.itp.new ./4zw9_331trp/4zw9_minimized_trp331_proa.itp
mv ./5eqi_333trp/5eqi_minimized_trp333_proa.itp.new ./5eqi_333trp/5eqi_minimized_trp333_proa.itp
mv ./5eqi_91asp/5eqi_minimized_asp91_proa.itp.new 5eqi_91asp/5eqi_minimized_asp91_proa.itp
mv ./4zw9_89asp/4zw9_minimized_asp89_proa.itp.new 4zw9_89asp/4zw9_minimized_asp89_proa.itp
```

-----------------------------------------------------------------------

```text
######## Wildtype
# To keep: 2034/2545 79.9%
# To skip: 511/2545 20.1%
############
# NEW: 
# To keep: 719/2461 29.2%
# To skip: 1742/2461 70.8%

python ../../code/comdyn.py --itp_file1 4zw9_wildtype/step6_proa.itp --itp_file2 5eqi_wildtype/5eqi_proa.itp

mv output.itp 4zw9_wildtype/4zw9_proa.itp

# To keep: 2034/2318 87.7%
# To skip: 284/2318 12.3%
python ../../code/comdyn.py --itp_file1 5eqi_wildtype/5eqi_proa.itp --itp_file2 4zw9_wildtype/4zw9_proa.itp
mv output.itp 5eqi_wildtype/5eqi_proa.itp
```

# 2. Do equilibration

First copy the calibration script:
```bash
for protein in 4zw9_331trp 5eqi_333trp 4zw9_89asp 5eqi_91asp 4zw9_wildtype 5eqi_wildtype
do
cp ../../code/run_calibration.sh ./${protein}/
cp ../../code/mdp/* ./${protein}/
done
```

Now run: 
```bash
current=`pwd`
source /usr/local/gromacs/bin/GMXRC
for protein in 4zw9_331trp 5eqi_333trp 78666 5eqi_91asp 4zw9_wildtype 5eqi_wildtype
do
cd ${protein}
bash ./run_calibration.sh
cd ${current}
done
```

# 3. Production run
Do this on Snellius:

```bash
current=`pwd`
for protein in 4zw9_331trp 5eqi_333trp 4zw9_89asp 5eqi_91asp 4zw9_wildtype 5eqi_wildtype
do
cd ${protein}
sbatch slurm_start_simulation8h_propermpi.sh
cd ${current}
done


current=`pwd`
for protein in 4zw9_331trp  4zw9_89asp  4zw9_wildtype
do
cd ${protein}
sbatch slurm_start_simulation8h_propermpi.sh
cd ${current}
done
```