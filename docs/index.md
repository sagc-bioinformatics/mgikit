---
title: Documenation and User Guide for MGIKIT
contributors: [Ziad Al-Bkhetan]
description: 
toc: false
---

## MGIKIT 
mgikit is a collection of tools used to demultiplex fastq files and generate demultiplexing and quality reports.

The toolkit includes the following commands:

### demultiplex
This command is used to demultiplex fastq files and assign the sequencing reads to their
associated samples. The tool requires the following mandatory input files to perform the
demultiplexing:
1. Fastq files (single/paired-end).
2. Sample sheet which contains sample indexes and their templates (will be explained in detail).

Simply, the tool reads the barcodes at the end of R2 (reveres) reads for paired-end reads input or the end of
R1 (forward) reads for single read input. Based on the barcode, it assigns the read to the relevant
sample allowing for mismatches less than a specific threshold. The tool outputs fastq files for each sample
as well as some summary reports that can be visualised through the MultiQC tool and mgikit plugin.

<hr/>

### template

This command is used to detect the location and form of the indexes within the read barcode. It simply goes though a small number of the reads and investiaget the number of matches with the indexes in th sammple sheet within each possible location in teh read barcode and consdering the indexes as is and their reverse complemntary. 

It reports matches for all possible combinations, and uses the read template that had the maximum number of matches. This process happens for each sample individually and therefore, the best matching template for each sample will be reported. 

Using this compehansive scan, teh tool can detect the templates for mixed libraries. 


### report

Under development!

<hr/>

## User Guide Table of Content

{% include section-navigation-tiles.html type="guides" %}

