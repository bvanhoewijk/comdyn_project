#!/bin/sh
addP=true
echo '#include "toppar/martini_v3.0.0.itp"'                     >gromacs/system.top
echo '#include "toppar/martini_v3.0.0_ions_v1.itp"'            >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_nucleobases_v1.itp"'     >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_phospholipids_v1.itp"'   >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_phospholipids_v1_matthieu.itp"' >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_small_molecules_v1.itp"' >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_solvents_v1.itp"'        >>gromacs/system.top
echo '#include "toppar/martini_v3.0.0_sugars_v1.itp"'          >>gromacs/system.top
echo '#include "toppar/martini_v3.0_sterols_v1.0.itp"'         >>gromacs/system.top
addP=false


for i in `ls gromacs/*.itp`;
do 
    ii=`basename $i`;
    echo "#include \"${ii}\"" >>gromacs/system.top;
done
echo "
[ system ]
; name
Martini system

[ molecules ]
; name        number" >>gromacs/system.top
for i in `grep ^ATOM step1_pdbreader.pdb |awk '{print substr($0, 73, 4)}'|uniq`
do
    chid=`grep ^ATOM step1_pdbreader.pdb|grep $i|head -1|awk '{print substr($0, 22, 1)}'`
    if [ $chid == " " ] || [ $addP = false ]
    then
        echo "${i} 1" >>gromacs/system.top
    else
        echo "${i}_${chid} 1" >>gromacs/system.top
    fi
done


cat step5_resname.str|uniq -c|awk '{print $2,$1}'>>gromacs/system.top
rm -rf step5_resname.str

