# Subsampling TH+ Amacrine data
This folder contains the scripts used for subsampling analysis of the TH+ amacrine sequencing data for Extended Data Figure 6d.  
A single representative experimental condition was chosen for the analysis.  

### Files you need
Make sure you have the .BAM file `TLS23A_Sorted.bam` in this directory with the README file.  
You can find the BAM alignment file in the posted repository for Light-Seq data.  

### Simulate Data.
The subsampling is done by creating a new .bam file using only a fraction of the reads, from 0.005 to 0.9.  
Each fraction condition is subsampled 5 times and the final count for each fraction is the mean +/- SD of the 5 simulations.  
Subsampling is done without replacement.  
First run:
```
python th_subsample.py
```
Which will generate a set of subsampled BAM files in the `simOut` folder. 

Next run:
```
python th_sortdata.py
python th_savedata.py
```
In that order, which will perform the same analysis of UMI deduplication on the simulated BAM files and sort them for the .csv file.  
