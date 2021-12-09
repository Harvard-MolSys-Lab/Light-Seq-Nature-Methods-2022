# Light-Seq Data Analysis

This repository contains code for analyzing Light-Seq sequencing data as presented in our [publication](TODO). Steps followed are outlined below. Steps were run either as batched jobs on the Harvard O2 server or locally on a 2018 Macbook Pro. 

## 1. Index separate and combined mouse and human genomes with STAR

[Human (Release 38)](https://www.gencodegenes.org/human/) and [mouse (Release M27)](https://www.gencodegenes.org/mouse/) assemblies and GFF3 files were downloaded. The mouse chromosomes were renamed to be `mchr` instead of `chr` so that the human and mouse genomes could be merged together for mapping the human/mouse cell mixing experiment data:

    $ sed -i 's/chr/mchr/g' GRCm39.primary_assembly.genome.fa
    $ sed -i -r 's/^chr/mchr/' gencode.vM27.annotation.gff3
    $ sed -i -r 's/ chr/mchr/' gencode.vM27.annotation.gff3

Then, merged assembly and GFF3 files were created:

    $ cat GRCh38.primary_assembly.genome.fa GRCm39.primary_assembly.genome.fa > GRCh38andGRCm39.primary_assembly.genome.fa
    $ cat gencode.v38.annotation.gff3 gencode.vM27.annotation.gff3 > gencode.v38andvM27.annotation.gff3

The [STAR aligner (Version 2.7.9a)](https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz) was downloaded from GitHub and used to index the human, mouse, and merged genomes:

    $ STAR-2.7.9a/source/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir 211101_GRCh38_Indexed --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v38.annotation.gff3 --sjdbOverhang 249
    $ STAR-2.7.9a/source/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir 211101_GRChm39_Indexed --genomeFastaFiles GRCm39.primary_assembly.genome.fa --sjdbGTFfile gencode.vM27.annotation.gff3 --sjdbOverhang 249
    $ STAR-2.7.9a/source/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir 211102_GRCh38andGRCm39_Indexed --genomeFastaFiles GRCh38andGRCm39.primary_assembly.genome.fa --sjdbGTFfile gencode.v38andvM27.annotation.gff3 --sjdbOverhang 249

## 2. Set up Conda environment

Next, you can set up your conda environment, as most scripts are written in Python. We recommend the following conda environment:

    $ conda create --name LSEnv python=3.7.5 pandas matplotlib seaborn PyTables Biopython scikit-image
    $ conda activate LSEnv
    $ conda install -c bioconda pysam umi_tools subread
    $ conda install -c https://conda.anaconda.org/biocore scikit-bio

Note: if samtools is not installed on your conda, you will need to install samtools (e.g. v1.9 or v1.12).  
With the LSEnv active:  

    $ conda install -c bioconda samtools==1.12

### (Optional) Set up the jupyter notebook kernel.

The following is optional if you want to set up the LSEnv as a jupyter notebook kernel.  
With LSEnv active:

    $ conda install ipykernel
    $ ipython kernel install --user --name=LSEnv

## 3. Set up directory structure

Indexed genomes should be placed in the main directory (where this README.md document is located), so that they can be accessed by the scripts CellMixing and RetinaTissue subfolder. If they are located elsewhere, you will need to update the scripts to call the correct locations. GFF3 annotation files should also be placed here.

Scripts for the cell mixing experiment are under the CellMixing subfolder, and scripts for the retina tissue experiment are under the RetinaTissue subfolder. [Raw imaging data](TODO) should be added into the inFiles subdirectory of the appropriate experiment. Scripts are written to have input fastq.gz files in the `inFiles` subdirectory and write output data to the `outFiles` subdirectory. You should start with an empty `outFiles` subdirectory, into which several output files for the experiment will be written. Next, the experiment-specific UMI extraction, genome mapping, transcript mapping, and analysis scripts can be run. 

## 4. Map and analyze cell mixing experiment

Once you have the indexed genomes and GFF3 files in the main directory and raw sequencing data files in the inFiles subdirectory, you can move to the CellMixing subdirectory to first extract UMIs:

    $ python3 ExtractTwoBarcodes.py

Next, UMI extracted reads can be mapped to the appropriate genomes:

    $ python3 MapToMergedAndSeparateHumanMouse.py

Note that for the cell mixing experiment, reads were mapped to both the combined human + mouse genome for discrimination analysis and to individual genomes for supplemental UMI estimates. Only reads that mapped uniquely a single time to the genome were considered. This means that all multi-mapped reads were excluded from further analysis and that the total number of barcoded molecules was likely substantially higher.

Now, sequences can be mapped to transcripts and deduped by the UMI + gene combo.

    $ python3 Dedup.py

Note that for mouse reads, genes that mapped to the ENSMUSG00000119584.1 and ENSMUSG00000064337.1 transcripts were excluded for further analysis, since they stalled the deduping.

GFP sequences can then be mapped:

    $ python3 AnalyzeGFPSequences.py

Once the sequencing files have been UMI extracted, Mapped, and UMI deduped. You should have several `_Dedup.bam` files in the `outFiles` folder. Several of the subsequent python scripts use the "pysam" package which requires an indexed '.bam' file. For convenience run:

    $ python3 samtoolsindexer.py

in the `outFiles` folder which will run the samtools on all the deduped files.  

Afterwards, the cellmixing numbers used for the publication can be analyzed by running:

    $ python3 cellmixing_1_parsegenes.py
    $ python3 cellmixing_2_calcTPM.py
    $ python3 cellmixing_3_plotTPM.py
    $ python3 cellmixing_4_parseUMIs.py
    
## 5. Map and analyze retina tissue experiment

Once you have the indexed genomes and GFF3 files in the main directory and raw sequencing data files in the inFiles subdirectory, you can also move to the RetinaTissue subdirectory to first extract UMIs:

    $ python3 ExtractThreeBarcodes.py

Next, UMI extracted reads can be mapped to the mouse genome:

    $ python3 MapToMouseUniqueOnly.py

Only reads that mapped uniquely a single time to the genome were considered. This means that all multi-mapped reads were excluded from further analysis and that the total number of barcoded molecules was likely substantially higher. 

Now, sequences can be mapped to transcripts and deduped by the UMI + gene combo.

    $ python3 Dedup.py

Note that for mouse reads, genes that mapped to the ENSMUSG00000119584.1 and ENSMUSG00000064337.1 transcripts were excluded for further analysis, since they stalled the deduping.

Now, tables of counts per gene for Light-Seq data can be generated (under the LightSeq subfolder):

    $ python3 VisualizeResultsLightSeq.py

Tables of simulated gene and cell counts per layer for Drop-Seq data can be generated (under the DropSeq subfolder):

    $ python3 VisualizeResultsDropSeq.py

Next R version 4.1.1 was used on RStudio to run [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) differential gene expression analysis to generate lists of significantly enriched (adjusted p value < 0.05) genes for each pair of layers for both Light-Seq and simulated Drop-Seq data (under their respective subfolders). These scripts were originally modeled off of R code written for the [Probe-Seq](https://elifesciences.org/articles/51452) method. For Light-Seq, a heatmap is also generated showing the genes for each layer that were significantly enriched compared to both other layers (shown in Figure 1 of the publication).

    > source("LightSeqDEA.R") 
    > source("LightSeqDEA.R") 

Finally, the plots shown in the main and supplementary figures of the publication can be generated:

    $ python3 PlotComparisons.py
