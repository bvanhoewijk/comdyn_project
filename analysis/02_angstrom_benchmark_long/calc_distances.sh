mkdir -p distances/5eqi
for item in 0.1nm 0.01nm 0.05nm
do
    for i in 1
    do
    tpr_file=5eqi/${item}/rep${i}/step7_production_rep${i}.tpr
    if [ ! -f ${tpr_file} ]; then
        echo "File ${tpr_file} does not exist. Continuing with next file."
        continue
    fi
    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM5_TM11_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_in plus com of group TM11_in'

    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM5_TM11_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_out plus com of group TM11_out'

    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM1_TM7_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_in plus com of group TM7_in'

    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM1_TM7_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_out plus com of group TM7_out'

    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM2_TM8_in_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_in plus com of group TM8_in'

    gmx distance \
        -s 5eqi/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 5eqi/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/5eqi/${item}_TM2_TM8_out_rep${i}.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_out plus com of group TM8_out'
    done
done


mkdir -p distances/4zw9
source /usr/local/gromacs/bin/GMXRC

for item in 0.1nm 0.01nm 0.05nm
do
    for i in 1
    do
    tpr_file=4zw9/${item}/rep${i}/step7_production_rep${i}.tpr
    if [ ! -f ${tpr_file} ]; then
        echo "File ${tpr_file} does not exist. Continuing with next file."
        continue
    fi
    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM5_TM11_in_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM5_in plus com of group TM11_in'

    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM5_TM11_out_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM5_out plus com of group TM11_out'

    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM1_TM7_in_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM1_in plus com of group TM7_in'

    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM1_TM7_out_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM1_out plus com of group TM7_out'

    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM2_TM8_in_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM2_in plus com of group TM8_in'

    gmx distance \
        -s 4zw9/${item}/rep${i}/step7_production_rep${i}.tpr \
        -f 4zw9/${item}/rep${i}/step7_production_rep${i}.xtc \
        -oall distances/4zw9/${item}_TM2_TM8_out_rep${i}.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM2_out plus com of group TM8_out'
    done
done

