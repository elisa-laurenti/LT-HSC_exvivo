### This is the github repository for the paper:
## Adaptation to ex vivo culture drives human haematopoietic stem cell loss of repopulation capacity in a cell cycle independent manner


### Fastq processing

``s1_process_fastq_files``

The code to process fastq files, aligning reads to the genome. Generation of count matrices.

### Quality control 

``s2_QC_and_data_preprocessing``

Filtering bad quality cells from the count matrices based on the quality of the sequencing of the individual cells.

### Seurat package

``s3_Seurat_analysis``

Analysis using the Seurat package for QA81, QA22 and mPB datasets.

### Differential gene expression analysis

``s4_DESeq2_analysis``

### Differential variable analysis

``s5_BASiCS_with_spikein``

Comparing the variance of genes between subsamples with the package BASiCS in R.

### Gene set variation analysis

``s6_GSVA_analysis``

Code for analysing pathway activation using the package GSVA.

### clustering along time course

``s7_DEGreport_analysis``

Code for using the package DEGreport to find patterns along a time course. 

### Measuring transcriptomic entropy

``s8_scEntropy_analysis``

Code for using the package scEntropy to measure transcriptome orderliness.

### supplementary analysis

``sup1_Scanpy_analysis``

The analysis of the package Scanpy on the data.

