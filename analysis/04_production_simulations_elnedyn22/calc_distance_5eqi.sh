#!/usr/bin/env bash
mkdir -p distances/5eqi
source /usr/local/gromacs/bin/GMXRC

for item in 5eqi 5eqi_asp91 5eqi_trp333
do
    for i in 1 2 3
    do
    tpr_file=${item}/rep${i}/step7_production_rep${i}.tpr
    xtc_file=${item}/rep${i}/step7_production_rep${i}.xtc
    if [ ! -f ${tpr_file} ]; then
        echo "File ${tpr_file} does not exist. Continuing with next file."
        continue
    fi
    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM5_TM11_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_in plus com of group TM11_in'

    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM5_TM11_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_out plus com of group TM11_out'

    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM1_TM7_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_in plus com of group TM7_in'

    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM1_TM7_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_out plus com of group TM7_out'

    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM2_TM8_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_in plus com of group TM8_in'

    gmx distance \
        -s ${tpr_file} \
        -f ${xtc_file} \
        -oall distances/5eqi/${item}_TM2_TM8_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_out plus com of group TM8_out'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM2_TM11_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM2_out plus com of group TM11_out'
    
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM2_TM11_in_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM2_in plus com of group TM11_in'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM5_TM8_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM5_out plus com of group TM8_out'
    
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM5_TM8_in_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM5_in plus com of group TM8_in'
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM8_TM9_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM8_out plus com of group TM9_out'
    
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM8_TM9_in_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM8_in plus com of group TM9_in'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM2_TM3_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM2_out plus com of group TM3_out'
    
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM2_TM3_in_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM2_in plus com of group TM3_in'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM1_TM1_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM1_in plus com of group TM1_out'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM5_TM5_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM5_in plus com of group TM5_out'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM7_TM7_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM7_in plus com of group TM7_out'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM11_TM11_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM11_in plus com of group TM11_out'

    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM2_TM2_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM2_in plus com of group TM2_out'
    # gmx distance \
    #     -s ${tpr_file} \
    #     -f ${xtc_file} \
    #     -oall distances/5eqi/${item}_TM8_TM8_in_out_rep${i}.com.xvg \
    #     -n config/5eqi_distance.ndx \
    #     -select 'com of group TM8_in plus com of group TM8_out'
    done
done
