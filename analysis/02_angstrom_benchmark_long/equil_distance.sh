%%bash
source /usr/local/gromacs/bin/GMXRC

# TM1 <-> TM7 out
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM1_TM7_out.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM1_out plus com of group TM7_out'


# TM1 <-> TM7 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM1_TM7_in.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM1_in plus com of group TM7_in'

# TM5 <-> TM11 out
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM5_TM11_out.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM5_out plus com of group TM11_out'

# TM5 <-> TM11 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM5_TM11_in.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM5_in plus com of group TM11_in'
# TM2 <-> TM8 out
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM2_TM8_out.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM2_out plus com of group TM8_out'

# TM2 <-> TM8 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/4zw9/step6.6_equilibration.tpr \
        -oall distances/4zw9/cg_equil_TM2_TM8_in.com.xvg \
        -n config/4zw9_distance.ndx \
        -select 'com of group TM2_in plus com of group TM8_in'

################################################################

# TM1 <-> TM7 out        
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM1_TM7_out.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_out plus com of group TM7_out'

# TM1 <-> TM7 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM1_TM7_in.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM1_in plus com of group TM7_in'

# TM5 <-> TM11 out
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM5_TM11_out.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_out plus com of group TM11_out'


# TM5 <-> TM11 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM5_TM11_in.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM5_in plus com of group TM11_in'


# TM2 <-> TM8 out
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM2_TM8_out.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_out plus com of group TM8_out'


# TM2 <-> TM8 in
gmx distance \
    -f /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.gro \
    -s /home/bas/projects/comdyn_project/analysis/04_production_simulations_elnedyn22/5eqi/step6.6_equilibration.tpr \
        -oall distances/5eqi/cg_equil_TM2_TM8_in.com.xvg \
        -n config/5eqi_distance.ndx \
        -select 'com of group TM2_in plus com of group TM8_in'
