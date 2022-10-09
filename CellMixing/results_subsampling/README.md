The subsampling script was performed on the deep sequencing run S436E.  
A BAM file from mapping to the merged human and mouse genome was first sorted then indexed:
```
samtools sort S436E_R1_MergedAligned.sortedByCoord.out.bam -o S436E_R1_MergedAligned.sortedByCoord.out_sorted.bam

samtools index S436E_R1_MergedAligned.sortedByCoord.out_sorted.bam
```
\
Before running  
```
cellmix_subsample.py
```
