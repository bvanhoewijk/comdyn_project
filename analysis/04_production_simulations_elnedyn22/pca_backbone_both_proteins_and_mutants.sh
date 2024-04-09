%%bash
source /usr/local/gromacs/bin/GMXRC

mkdir -p pca/out


#### Wildtype
gmx trjcat -f \
./4zw9/rep1/step7_production_rep1.xtc \
./4zw9/rep2/step7_production_rep2.xtc \
./4zw9/rep3/step7_production_rep3.xtc \
-o pca/4zw9_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/4zw9_everything_all_rep.xtc \
    -s 4zw9/rep1/step7_production_rep1.tpr \
    -n 4zw9/index.bb.ndx \
    -o pca/4zw9_eigenval.everything.bb.xvg \
    -v pca/4zw9_eigenvec.everything.bb.trr


#### Mutant 1
gmx trjcat -f \
./4zw9_asp89/rep1/step7_production_rep1.xtc \
./4zw9_asp89/rep2/step7_production_rep2.xtc \
./4zw9_asp89/rep3/step7_production_rep3.xtc \
-o pca/4zw9_asp89_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/4zw9_asp89_everything_all_rep.xtc \
    -s 4zw9_asp89/rep1/step7_production_rep1.tpr \
    -n 4zw9_asp89/index.bb.ndx \
    -o pca/4zw9_asp89_eigenval.everything.bb.xvg \
    -v pca/4zw9_asp89_eigenvec.everything.bb.trr

#### Mutant 2
gmx trjcat -f \
./4zw9_trp331/rep1/step7_production_rep1.xtc \
./4zw9_trp331/rep2/step7_production_rep2.xtc \
./4zw9_trp331/rep3/step7_production_rep3.xtc \
-o pca/4zw9_trp331_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/4zw9_trp331_everything_all_rep.xtc \
    -s 4zw9_trp331/rep1/step7_production_rep1.tpr \
    -n 4zw9_trp331/index.bb.ndx \
    -o pca/4zw9_trp331_eigenval.everything.bb.xvg \
    -v pca/4zw9_trp331_eigenvec.everything.bb.trr

#### Wildtype
gmx trjcat -f \
./5eqi/rep1/step7_production_rep1.xtc \
./5eqi/rep2/step7_production_rep2.xtc \
./5eqi/rep3/step7_production_rep3.xtc \
-o pca/5eqi_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_everything_all_rep.xtc \
    -s 5eqi/rep1/step7_production_rep1.tpr \
    -n 5eqi/index.bb.ndx \
    -o pca/5eqi_eigenval.everything.bb.xvg \
    -v pca/5eqi_eigenvec.everything.bb.trr


#### Mutant 1
gmx trjcat -f \
./5eqi_asp91/rep1/step7_production_rep1.xtc \
./5eqi_asp91/rep2/step7_production_rep2.xtc \
./5eqi_asp91/rep3/step7_production_rep3.xtc \
-o pca/5eqi_asp91_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_asp91_everything_all_rep.xtc \
    -s 5eqi_asp91/rep1/step7_production_rep1.tpr \
    -n 5eqi_asp91/index.bb.ndx \
    -o pca/5eqi_asp91_eigenval.everything.bb.xvg \
    -v pca/5eqi_asp91_eigenvec.everything.bb.trr

#### Mutant 2
gmx trjcat -f \
./5eqi_trp333/rep1/step7_production_rep1.xtc \
./5eqi_trp333/rep2/step7_production_rep2.xtc \
./5eqi_trp333/rep3/step7_production_rep3.xtc \
-o pca/5eqi_trp333_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_trp333_everything_all_rep.xtc \
    -s 5eqi_trp333/rep1/step7_production_rep1.tpr \
    -n 5eqi_trp333/index.bb.ndx \
    -o pca/5eqi_trp333_eigenval.everything.bb.xvg \
    -v pca/5eqi_trp333_eigenvec.everything.bb.trr


# Wildtype
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 4zw9/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/4zw9_eigenval.everything.bb.xvg \
    -v pca/4zw9_eigenvec.everything.bb.trr \
    -f ./4zw9/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/4zw9_rep${i}.everything.bb.xvg -first 1 -last 2 -n 4zw9/index.bb.ndx 2> log.txt
done

# Mutant 1
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 4zw9_asp89/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/4zw9_asp89_eigenval.everything.bb.xvg \
    -v pca/4zw9_asp89_eigenvec.everything.bb.trr \
    -f ./4zw9_asp89/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/4zw9_asp89_rep${i}.everything.bb.xvg -first 1 -last 2 -n 4zw9_asp89/index.bb.ndx 2> log.txt
done

# Mutant 2
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 4zw9_trp331/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/4zw9_trp331_eigenval.everything.bb.xvg \
    -v pca/4zw9_trp331_eigenvec.everything.bb.trr \
    -f ./4zw9_trp331/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/4zw9_trp331_rep${i}.everything.bb.xvg -first 1 -last 2 -n 4zw9_trp331/index.bb.ndx 2> log.txt
done


#### Wildtype
gmx trjcat -f \
./5eqi/rep1/step7_production_rep1.xtc \
./5eqi/rep2/step7_production_rep2.xtc \
./5eqi/rep3/step7_production_rep3.xtc \
-o pca/5eqi_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_everything_all_rep.xtc \
    -s 5eqi/rep1/step7_production_rep1.tpr \
    -n 5eqi/index.bb.ndx \
    -o pca/5eqi_eigenval.everything.bb.xvg \
    -v pca/5eqi_eigenvec.everything.bb.trr \
    -ref


#### Mutant 1
gmx trjcat -f \
./5eqi_asp91/rep1/step7_production_rep1.xtc \
./5eqi_asp91/rep2/step7_production_rep2.xtc \
./5eqi_asp91/rep3/step7_production_rep3.xtc \
-o pca/5eqi_asp91_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_asp91_everything_all_rep.xtc \
    -s 5eqi_asp91/rep1/step7_production_rep1.tpr \
    -n 5eqi_asp91/index.bb.ndx \
    -o pca/5eqi_asp91_eigenval.everything.bb.xvg \
    -v pca/5eqi_asp91_eigenvec.everything.bb.trr

#### Mutant 2
gmx trjcat -f \
./5eqi_trp333/rep1/step7_production_rep1.xtc \
./5eqi_trp333/rep2/step7_production_rep2.xtc \
./5eqi_trp333/rep3/step7_production_rep3.xtc \
-o pca/5eqi_trp333_everything_all_rep.xtc

printf "BB BB\n" | gmx covar \
    -f pca/5eqi_trp333_everything_all_rep.xtc \
    -s 5eqi_trp333/rep1/step7_production_rep1.tpr \
    -n 5eqi_trp333/index.bb.ndx \
    -o pca/5eqi_trp333_eigenval.everything.bb.xvg \
    -v pca/5eqi_trp333_eigenvec.everything.bb.trr


# Wildtype
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 5eqi/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/5eqi_eigenval.everything.bb.xvg \
    -v pca/5eqi_eigenvec.everything.bb.trr \
    -f ./5eqi/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/5eqi_rep${i}.everything.bb.xvg -first 1 -last 2 -n 5eqi/index.bb.ndx 2> log.txt
done

# Mutant 1
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 5eqi_asp91/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/5eqi_asp91_eigenval.everything.bb.xvg \
    -v pca/5eqi_asp91_eigenvec.everything.bb.trr \
    -f ./5eqi_asp91/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/5eqi_asp91_rep${i}.everything.bb.xvg -first 1 -last 2 -n 5eqi_asp91/index.bb.ndx 2> log.txt
done

# Mutant 2
for i in 1 2 3
do
echo "Replicate ${i}"
printf "BB BB\n" | gmx anaeig \
    -s 5eqi_trp333/rep${i}/step7_production_rep${i}.tpr  \
    -eig pca/5eqi_trp333_eigenval.everything.bb.xvg \
    -v pca/5eqi_trp333_eigenvec.everything.bb.trr \
    -f ./5eqi_trp333/rep${i}/step7_production_rep${i}.xtc \
    -2d pca/out/5eqi_trp333_rep${i}.everything.bb.xvg -first 1 -last 2 -n 5eqi_trp333/index.bb.ndx 2> log.txt
done