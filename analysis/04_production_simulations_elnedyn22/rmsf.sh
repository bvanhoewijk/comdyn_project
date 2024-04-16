#!/usr/bin/env bash
source /usr/local/gromacs/bin/GMXRC

# GROUP 16 is the Backbone
####################### RMSD 4ZW9
for item in 4zw9 4zw9_asp89 4zw9_trp331
do
    # RSMD
    for i in 1 2 3
    do
    tpr_file=${item}/rep${i}/step7_production_rep${i}.tpr
    xtc_file=${item}/rep${i}/step7_production_rep${i}.xtc
    if [ ! -f ${tpr_file} ]; then
        echo "File ${tpr_file} does not exist. Continuing with next file."
        continue
    fi
    printf "16" | gmx rmsf \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -n ${item}/index.bb.ndx \
        -o rmsf/${item}_rep${i}.xvg \
        -b 1000 \
        -fit
    done
done

####################### RMSD 5EQI
for item in 5eqi 5eqi_asp91 5eqi_trp333
do
    # RSMD
    for i in 1 2 3
    do
    tpr_file=${item}/rep${i}/step7_production_rep${i}.tpr
    xtc_file=${item}/rep${i}/step7_production_rep${i}.xtc
    if [ ! -f ${tpr_file} ]; then
        echo "File ${tpr_file} does not exist. Continuing with next file."
        continue
    fi
    printf "16" | gmx rmsf \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -n ${item}/index.bb.ndx \
        -o rmsf/${item}_rep${i}.xvg \
        -b 1000 \
        -fit
    done
done
