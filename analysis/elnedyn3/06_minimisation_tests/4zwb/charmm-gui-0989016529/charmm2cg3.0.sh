#!/bin/sh
mkdir -p gromacs/toppar
i=$1
molname=$2
charge=$3
elastic=$4

if [[ $charge == "neutral" ]]
then
    charge="-nt"
else
    charge=" "
fi

if [[ $elastic == "elastic" ]]
then
    elastic="-elastic"
else
    elastic=" "
fi

rm *.log
echo "./martinize2 -f ${i}.charmm.pdb -o ${i}.top -x ${i}.cg.pdb -ff martini3001 $charge $elastic -dssp /home/charmm-gui/local/dssp-3.1.4/mkdssp -p backbone -maxwarn 1 -mutate HSD:HIS -mutate HSP:HIH" > ${i}.log
./martinize2 -f ${i}.charmm.pdb -o ${i}.top -x ${i}.cg.pdb -ff martini3001 $charge $elastic -dssp /home/charmm-gui/local/dssp-3.1.4/mkdssp -p backbone -maxwarn 1 -mutate HSD:HIS -mutate HSP:HIH >> ${i}.log 2>&1
mv molecule_0.itp ${i}.itp

#f=${i/.charmm.pdb/.cg.pdb}
f=`echo ${i}.charmm.pdb|sed s/.charmm.pdb//`
segid=`cat -v ${i}.charmm.pdb|grep ^ATOM|awk '{print substr($0,73,4)}'|head -1`
echo $f $segid
 
awk -v ch="$segid" '{if($1=="ATOM") printf("%s      %s\n",substr($0, 0, 66),ch);else print}' $f.cg.pdb>temp
mv -f temp $f.cg.pdb
#/usr/local/opt/gnu-sed/libexec/gnubin/sed -i s/molecule_0/$molname/ $f.itp
#sed -i -e "s/molecule_0/$molname/" -e "s/1000[[:space:]]1000[[:space:]]1000/POSRES_FC\ POSRES_FC\ POSRES_FC/" -e "s/#ifdef POSRES/#ifdef POSRES\n#ifndef POSRES_FC\n#define POSRES_FC 1000.00\n#endif/" $f.itp
sed -i "s/molecule_0/$molname/" $f.itp
/home/charmm-gui/local/bin/python3_conda /home/charmm-gui/www/lib/normang.py $f.itp
mv $f.itp gromacs
