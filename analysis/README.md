# Folder setup:

# Elnedyn3: Benchmarking using the Martini3 forcefield.

- Note: 4ZW9 protein breaks during the simulations. 
- For the thesis the data from Elnedyn22 was used.

# Elnedyn22: Benchmarking and produciton simulations using the Martini v2.2 forcefields


## 01_angstrom_benchmark_short
Description: Benchmark analysis using seven thresholds for removal of rubber bands of 2microseconds

Analysis notebooks:

- 01_distance_analysis.ipynb
  - Inter helix distances
  - KDE plots

- 02_pca_same_space.ipynb
  - PCA analysis of general movement of the protein 

- 03_rmsd.ipynb
  - RMSD calculations
 

- 04_radius_of_gyration.ipynb
  - Radius of gyration calculations
  


## 02_angstrom_benchmark_long
Description: Benchmark analysis using 4 thresholds for removal of rubber bands of 50 microseconds

- 01_distance_analysis_longbenchmark.ipynb
  - inter TM distnaces over time
  - Also adds starting crystal location

- inital_distances.tsv
  - Initial distances between TM helics as calculated by `init_distance2tsv.py`

## 03_make_mutants
Description: Convenience script for making protein structure mutations

## 04_production_simulations

- 01_distance_analysis.ipynb
- 02_pca_same_space.ipynb
- 03_rmsd.ipynb
- 04_radius_of_gyration.ipynb
- 05_rmsf.ipynb
- 06_distance_overlap_shift.ipynb
- 07_pca_overlap.ipynb

- *.tgz: Wildtype/Mutant systems. (uncalibrated)

