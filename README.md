---

![SAGC-Bioinformatics](docs/assets/SAGC-logo-hover.png)

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

This command is used to detect the location and form of the indexes within the read barcode. It simply goes through a small number of the reads and investigates the number of matches with the indexes in the sample sheet within each possible location in the read barcode and considering the indexes as is and their reverse complementary. 

It reports matches for all possible combinations and uses the read template that had the maximum number of matches. This process happens for each sample individually and therefore, the best matching template for each sample will be reported. 

Using this comprehensive scan, the tool can detect the templates for mixed libraries. 


### report

This command is to merge demultiplexing and quality reports from multiple lanes into one comprehensive report for MultQC reports visualisation.

<hr/>

## Installation

You can use the static binary under bins directly, however, if you like to build it from the source code:

You need to have Rust and cargo installed first, check rust [documenation](https://doc.rust-lang.org/cargo/getting-started/installation.html)

```bash
git clone https://github.com/sagc-bioinformatics/mgikit.git
cd mgikit
cargo build --release
```



## User Guide

Please checkout the [documeantion](https://sagc-bioinformatics.github.io/mgikit/)


## Commerical Use

Please contact us if you want to use the software for commercial purposes.
