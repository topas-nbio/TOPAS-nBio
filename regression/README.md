### A regression system for TOPAS-nBio-dev                          ###
### Jose Ramos-Mendez and Naoki D. Kondo, Thongchai A. M. Masilela    #
### Department of Radiation Oncology                                  #
### University of California San Francisco.                           #
### Jose.RamosMendez@ucsf.edu                                         #
### Mar 2024.                                                         #
### ############################################################### ###

#### Aim
Perform a regression test for TOPAS-nBio based on examples with existing published reference data.
#### Test
1. DBSCAN. Quantifies the ratio between SSB to DSB for protons using the clustering algorithm DBSCAN. Reference data from previous Geant4-DNA version (Francis et al., 2017) and digitized data from PARTRAC and measured is provided.

2. GvalueStepByStep. Quantifies the G-value as a function of the time for fast electrons of 1 MeV. Reference data from Wang et al., 2018 at the shortest times (7 ps) is provided.

3. LET. Quantifies restricted LET. Reference data from PARTRAC is provided.

4. NanodosimetryI. Quantifies the conditional ionization cluster size distribution (nu > 0) from carbon ions of 88 MeV incident in a single cylinder. Measured data obtained in gas from Hilgers et al., 2017 is provided.

5. NanodosimetryII. Quantifies F2 vs M1, where F2 is the cumulative probability of having ionization clusters of size 2 or bigger, and M1 is the fist moment of the unconditional ionization cluster size distribution. Measured data obtained from Conte et al. 2017 is provided.

6. NanodosimetryIII. Quantifies ionization cluster size distribution produced in randomly oriented cylinders for ions from proton to oxygen ions. Calculated reference data from Ramos-Mendez et al., 2017 is provided.

7. GValueIRT. Quantifies the G-value as a function of the time for fast electrons of 1 MeV. 

8. FrickeIRT. Quantifies the G-value of Fe^3+ which comes from the oxidation of Fe^2+.
This example must give a value of around 15.5 +- 0.1 reported by the ICRU.

9. GValue_LET-IRT. LET-dependent G values for e-, p and alpha at selected energies

10. GValue_LET-SBS. LET-dependent G values for e-, p and alpha at selected energies

11. GvalueIRT-Temperature. Temperature-dependent G values for fast electrons at T < 200C

12. GvalueIRT_H2O2. OH-scavenger-dependent G value for H2O2 within microsecond time range

13. GvalueIRT_H. H-scavenger-dependent G value for H within microsecond time range

14. DNASSBPulsed. Quantifies the number of SSBs induced on a DNA plasmid as a function of DMSO concentration for a pulsed beam of 250 keV electrons

#### Use
Each directory has TOPAS parameters with the follow pattern:
- `mainABCD.txt`, where ABCD is Topas, Opt2, Opt4 or Opt6 that refer to used physics list.
- `depFileN.txt`, with N=1,... dependence file used by the main file.
- `inputfiles.txt`. Lists the main file (s) to be submitted as a job.
- `tcsh script files`. Each directory has three tcsh files to submit the example to e.g. a cluster system or locally. For instance, to submit 5 simulation jobs in the local system (different random seeds) use

```bash
tcsh submitLocally.sh 5
```
A directory named `run` is created which contains a sub-directory, its name contains the current date. Later, compare between simulation results with the python script in `analysis/analysis.py`

```bash
python analysis/analysis.py run/2020July/mainTopas run/2020July/mainOpt2 --sut_label Topas --ref_label Opt2
```
Look for the images in the directory `results/`. A table which contains averaged execution time per CPU is also available.

#### Optional
Run the Python script inside the folder Summary/tex_openTOPAS. This will copy and paste all the regression test images into a new directory called Summary/openTOPAS. Run the .tex file in order to generate a PDF summarising all the results.


