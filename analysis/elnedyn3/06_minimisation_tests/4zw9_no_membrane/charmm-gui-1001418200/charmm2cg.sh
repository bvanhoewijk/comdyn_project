#!/bin/sh
mkdir  -p gromacs/toppar
ff=$1
i=$2
charge=$3
awk '{if($0~"^ATOM") {printf("%s%5d%s\n",substr($0,1,6),a+1,substr($0,12,90));a++} else print}' $i >temp.pdb
mv temp.pdb $i
f=`echo $i|sed s/.charmm.pdb//`
segid=`cat -v $i|grep ^ATOM|awk '{print substr($0,73,4)}'|head -1`
echo $f $segid
./mkdssp -i $f.charmm.pdb -o $f.charmm.dssp

if [[ $ff == "elnedyn" ]]
then
    ff="elnedyn22"
elif [[ $ff == "elnedynp" ]]
then
    ff="elnedyn22p"
fi

if [[ $charge == "neutral" ]]
then
    echo "neutral"
    python martinize.py -f $f.charmm.pdb -o $f.cg.top -x $f.cg.pdb -ss $f.charmm.dssp -ff $ff -p Backbone -name  $segid -cys auto -nt
else
    echo "charged"
    python martinize.py -f $f.charmm.pdb -o $f.cg.top -x $f.cg.pdb -ss $f.charmm.dssp -ff $ff -p Backbone -name  $segid -cys auto 
fi
 
awk -v ch="$segid" '{if($1=="ATOM") printf("%s      %s\n",$0,ch);else print}' $f.cg.pdb>temp
mv -f temp $f.cg.pdb

if [[ $ff == "martini22p" || $ff == "elnedyn22p" ]]
then
    echo "p"
    #in martinize.py martini22p output, the dummy atoms are all named SCD
    #ATOM      4  SC1 ASN A   2      35.659 -21.603  13.346  1.00  0.00      PROA
    #ATOM      5  SCD ASN A   2      35.551 -21.515  13.363  1.00  0.00      PROA
    #ATOM      6  SCD ASN A   2      35.767 -21.690  13.328  1.00  0.00      PROA
    #this awk script renames SCD to SCP, SCN ... based on the previous sidechain atom
    awk '{
        if(substr($0,0,4) == "ATOM")
        {
            atomname = substr($0,13,4)
            resname  = substr($0,17,4)
            if(atomname != " SCD")
            {
                print
                atomprev = atomname
                next
            }
            if(atomprev == " SCP")
                atomname = " SCN"
            else
                atomname = " SCP"
    
            if((resname == " ASP"||resname == " GLU") && atomname == " SCP")
                atomname = " SCN"
            atomprev = atomname
            printf("%s%s%s%s\n",substr($0,0,12), atomname, resname, substr($0,21,100))
        } else {
            print
        }
    
    }' $f.cg.pdb >temp
    mv -f temp $f.cg.pdb
fi

mv *.itp gromacs

