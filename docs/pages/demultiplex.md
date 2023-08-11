---
title: Instructions for MGIKIT demultiplexing
contributors: [Ziad Al-Bkhetan]
description: 
toc: true
type: guides
---

### Demultiplexing functionality
This command is used to demultiplex fastq files and assign the sequencing reads to their
associated samples. The tool requires the following mandatory input files to perform the
demultiplexing:
1. Fastq files (single/paired-end).
2. Sample sheet which contains sample indexes and their templates (will be explained in detail).

Simply, the tool reads the barcodes at the end of R2 (reveres) reads for paired-end reads input or the end of
R1 (forward) reads for single read input. Based on the barcode, it assigns the read to the relevant
sample allowing for mismatches less than a specific threshold. The tool outputs fastq files for each sample
as well as some summary reports that can be visualised through the MultiQC tool and mgikit plugin.

### Mandatory input files

**Fastq input file**

+ **`-f or --read1`**: the path to the forward reads fastq file for both paired-end and single-end input data.

+ **`-r or --read2`**: the path to the reverse reads fastq file.

+ **`-i or --input`**: the path to the directory that contains the input fastq files. 

Either `-i` or `-f/-r`, `-f` should be provided for a run.

**Input sample sheet**

+ **`-s or --sample-sheet`**: the path to the sample sheet file. 

More details are available below on the sample sheet format and preparation.

### Other Parameters

+ **`-o or --output`**: The path the output directory.

The tool will create the directory if it does not exist
or overwrite the content if the directory exists and the parameter `--force` is used. The tool will exit
with an error if the directory exists, and `--force` is not used. If this parameter is not provided, the tools
will create a directory (in the working directory) with a name based on the date and time
of the run as follows `mgiKit_Y-m-dTHMS`. where `Y`, `m`, `d`, `H`, `M`, and `S` are the date and time format.

+ **`--reports`**: The path of the output reports directory.

By default, the tool writes the files of the run reports in the same output directory as the
demultiplexed fastq files (`-o` or `--output` parameter). This parameter is used to write the reports in
a different folder as specified with this parameter.

+ **`-m or --mismatches`**: The default value is 1. The number of mismatches allowed when
matching reads’ barcode with sample indexes. 

This number should be less than the minimal Hamming
distance between any barcodes of two samples. In the case of dual index demultiplexing (i7 and i5), a
read will be assigned to a sample if the sum of the mismatches between the read barcode and both
indexes is less than or equal to the value of this parameter. In the case of a single index, the
mismatches with the single index should be less or equal to this parameter.

+ **`--disable-illumina`**: Output reads' header in MGI format.

This option is to disable the default behaviour of the tool that outputs read files using Illumine format (for read headers and file naming). More details are below.

+ **`--keep-barcode`**: keep the barcode at the end of the demultiplexed read.

By default, the tool trims the barcode sequence at the end of the read sequence. This can be disabled using this flag and the demultiplexed reads will contain the barcode at the tail of the read2 for paired-read sequencing or the tail of read1 for single-read sequencing.

+ **`--template`**: The general template of the index locations in the read barcode. details are in the sample sheet preparation. if all samples in the sample sheet use the same template, the general template can be passed as a parameter instead of having it in the sample sheet for each sample. The general template will be used for all samples with the combination of `--i7-rc` and `--i5-rc`.

+ **`--i7-rc`**: Use the reverse complementary form of i7 when matching with the barcode. 

This option should be used together with the general template and will be applied to all samples.

+ **`--i5-rc`**: Use the reverse complementary form of i5 when matching with the barcode. 

This option should be used together with the general template and will be applied to all samples.

+ **`--lane`**: Lane number such as `L01`.

This parameter is used to provide the lane number when the parameter `-i` or `--input` is not
provided. The lane number is used for QC reports and it is mandatory when Illumina format is
requested for file naming.

+ **`--instrument`**: The id of the sequncing machine. 

This parameter is used to provide the instrument id when the parameter `-i` or `--input`
is not provided. The parameter is mandatory when Illumina format is requested for read header and
file naming.

+ **`--run`**: The run id. It is taken from Bioinf.csv as the date and time of starting the run.

This parameter is used to provide the run id when the parameter `-i` or `--input` is not provided. The parameter is mandatory when Illumina format is requested for read header and file naming.

+ **`--merged-reads`**: The default value is `10,000`. The number of reads that will be concatenated in one string and then written at once with one write command.

+ **`--writing-buffer`**: The default value is `1000`. The number of merged reads in the buffer which will be written to the output file.

+ **`--comprehensive-scan`**: Enable comperhansive scan. 

This parameter is only needed when having a mixed library dataset (different locations for the indexes in the read barcode for some samples).
The default behaviour of the tool is to stop comparing the barcodes with
samples’ indexes when it finds a match. This flag will force the tool to keep comparing with all other
samples to make sure that the read matches with only one sample. In a normal scenario, the read
should match with only one sample, however, there is a chance that the read matches with multiple
samples if the allowed number of mismatches is greater than the minimum hamming distance
between the indexes, or the samples have different templates.

+ **`--undetermined-label`**: The default value is `Undetermined`. The label of the file that contains the
undermined reads which could not be assigned to any samples.

+ **`--ambiguous-label`**: The default value is `ambiguous`. The label of the file that contains the ambiguous reads.
The ambiguous reads are the reads that can be assigned to multiple samples. This can happen when
the number of allowed mismatches is high.

+ **`--force`**: this flag is to force the run and overwrite the existing output directory if exists.

### Understanding input files
MGI sequencing machine output a directory for the run (flowcell_id) with a subdirectory for each lane (L01, L02 ..) depending on the machine.
The input fastq files can be provided to the tool in two ways:

1.  Using `-f` and `-r` parameters which will be referring to the path to `R1` and `R2` respectively for paired-end or `-f` for single end fastq.

2. Using `-i` or `--input` parameter which refers to the path to the lane subdirectory in the sequencing output directory (or the directory that contains the fastq files if the data is obtained from somewhere else). In this case, the tool will search for the file that ends with `_read_1.fq.gz` and `_read_2.fq.gz` as forward and reverse reads respectively and if no reverse read file is found, the tool considers the run as a single end run. 

### Output fastq format

**File naming**

1. Illumina format (default format)

    The fils will have the following pattern `SAMPLEID_S{1-n}_L0{1,2,3,4}_R{1,2}_001.fastq.gz`.

    For example, `21-10233_S1_L01_R1_001.fastq.gz` and `21-10233_S1_L01_R2_001.fastq.gz`.
    <div align="center">
    
    ![read-headers-figure](/docs/assets/file-naming.png)
    
    </div>
    
    Where:
        
    1. Sample ULN.
    
    2. Sample order in the sample sheet.
    
    3. Lane number.
    
    4. R1 and R2 for paired end sequencing, or R1 for single end sequencing.
    
    5. The remaining part is fixed for all files.

2. MGI format 

The files will have the following pattern `SAMPLEID_L0{1,2,3,4}_R{1,2}.fastq.gz`.

For example, `21-10233_L01_R1.fastq.gz` and `21-10233_L01_R2.fastq.gz`.



The lane number is an input parameter and `R1` and `R2` are for forward and reverse read.

**Read header**

1. Illumina format (default format)

2. MGI format

Please see the details of both headers and the conversion in the figure below:

<div align="center">

![read-headers-figure](docs/assets/read-header.png)

</div>

Illumina formatting requires Lane number, instrument id and run id. The three requirements can be provided using their parameters `--lane`, `--instrument`, and `--run` as described above. However, if the `--input` was used, the tools will look for `BioInfo.csv` file that is generated by MGI sequencers in the input directory and extract the instrument id and run id from it. The run id will be the date and time of the run start ("YMDHmS" format). It will also check the second element in the name of the input fastq file (after splitting it by `_`) if it contains the substring `L0*`, if yes, the lane number will be taken from the fastq file name.

`--lane`, `--instrument`, and `--run` will be prioritised over the information in the `BioInfo.csv` file if both were provided.

### Sample sheet format and preparation.

For the tool to perform demultiplexing, it needs to know the indexes of each sample to match them with the barcodes at the end of the read sequence as well as where to look for each index in the barcode. We refer to the location of the indexes within the barcode by the barcode template. For example
This information can be provided to the tool in two different ways:
1. Sample sheet with all information about the indexes and their template.
This is the general way of passing this information to the tool and it works for all scenarios.
2. Sample sheet with information about sample indexes and a general template passed via a separate parameter.
This approach can be used if all samples have the same index template, which is the common scenario.



The sample sheet is a tab-delimited file that contains sample information. 

This file may contain all or some of these columns:

+ **`sample_id`** (**Mandatory**): a unique identifier that will be used to refer to the specific sample.

+ **`i7`** (**Mandatory**): The nucleotide sequence for the i7 index for the associated sample. This will be used to demultiplex the reads by comparing it to the index found in the read barcode.

+ **`i5`** (**Optional**): The nucleotide sequence for the i5 index for the associated sample. this will be used to demultiplex the reads by comparing it to the index found in the read barcode when using dual indexes.

+ **`template`** (**Optional**): This column should contain the template of the barcode for the specific sample. This allows doing demultiplexing for samples from different libraries where the templates are different. If all samples have the same template, this column can be ignored, and a general template should be passed in a separate parameter. See `--template` parameter. More details are below.

+ **`i7_rc`** (**Optional**): Takes values from 0 or 1. If the value is 0, the i7 in the sample sheet will be used as is to compare with the read barcode. If the value is 1, the tool will compare the index found in the read barcode to the reverse complementary of i7 in the sample sheet. If the template was not provided in the sample sheet (general template is used), this parameter will be ignored, and the user has to provide this parameter (`--i7-rc`) separately.

+ **`i5_rc`** (**Optional**): The same as i7_rc but applies to the i5 index.

+ **`job_number`** (**Optional**): It is an id to group the samples that are from the same project for the cases when a run contains samples from multiple projects. The demultiplexer will generate demultiplexing and quality reports for each project and the whole run. It can be ignored if the run has samples for the same project or if the project-based reports are not needed.


**Barcode template**

To understand how to use the demultiplexing tool, it is important to understand the structure of the input data and how to provide the correct parameters for the analysis.

The sequenced reads obtained by the MGI sequencing machine contain a string of nucleotides at the tail of read2 for paired-end sequencing or the tail of read1 for single-end sequencing. This substring is referred to as the read barcode which contains the indexes of the samples, single (i7) or dual (i7 and i5) indexes. It also includes the Unique Molecular Identifier (UMI) in some cases. 

The demultiplexer tool looks at this read barcode and tries to match the indexes to a subsequence within the barcode to assign the read to a specific individual. In order to accomplish that, the tool needs to know where to look at the barcode to match with the index and from where to extract the UMI. This information is provided to the tool through the template parameter.
The template parameter is a combination of four possible components of (i7*, i5*, um*, --*) separated by a colon `:`. Where:

+ `i7`: the region where to expect the index i7 within the barcode.
+ `i5`: the region where to expect the index i5 within the barcode
+ `um`: the region where to expect the UMI within the barcode
+ `--`: a discarded region that is not used.
+ `*` : should be replaced by a number representing the length of the relevant index or UMI.


*Examples of templates*

**Example 1:** The template i58:um8:i78 means:

1. The length of the barcode at the end of read2 (paired-end) or at the end of read1 (single-end) is 24 bp (8 + 8 + 8).

2. The last 8 bp of this barcode contains i7.

3. The middle 8 bp contains the UMI.

4. The first 8 bp contains i5.

<div align="center">

![template-example-1](docs/assets/template-example-1.png)

</div>

**Example 2:** The template --2:i58:--2:i78 means:

1. The length of the barcode at the end of read2 (paired-end) or at the end of read1 (single-end) is 20 bp (2 + 8 + 2 + 8).

2. The last 8 bp of this barcode contains i7.

3. The two base pairs before the i5’s 8 bp are not used and should be ignored.

4. The 8 base pairs before, contain the i5

5. The first 2 bp before should be ignored as well.

<div align="center">

![template-example-2](docs/assets/template-example-2.png)

</div>

Note that in this example the direction of the arrows is to show that these indexes are reverse complementary in the reads, therefore, this should be accounted for when demultiplexing using other parameters as explained in this documentation later.

**Example 3:** The template --2:i78 means:

1. The length of the barcode at the end of read2 (paired-end) or at the end of read1 (single-end) is 10 bp (2 + 8).

2. The last 8 bp of this barcode contains i7.

3. The first two base pairs before the i7’s 8 bp are not used and should be ignored.

4. There is no i5 in this template as it is a single index run.


After demultiplexing, the barcode will be trimmed by default including all parts mentioned in the template. This can be disabled using the parameter `--keep-barcode`.

Templates and indexes forms can be provided by the user, however, the command `template` can detect the barcode template and the form of the indexes for the run. 


### Reports

The demultiplex command generates multiple reports with file names that start with the flowcell and lane being demultiplexed:

1. `flowcell.L0*.mgikit.info`

This report contains the number of reads per sample respectively to each possible mismatch. It has (2 + allowed mismatches during demultiplexing) columns. 
For example:

<div align="center">

| **sample** | **0-mismatches** | **1-mismatches** |
|:------------:|:------------------:|:------------------:|
|     S01    |     3404         |      5655        |

</div>

This means that there was only one mismatch allowed during this execution and the sample S01 has 3404 reads with indexes matching perfectly and 5655 reads with indexes that differ by 1 base compared to the indexes provided in the sample sheet.

This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

2. `flowcell.L0*.mgikit.general`

This file contains summary information related to the cluster count and quality scores, summarised for each sample as well as at the whole lane scale. This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

3. `flowcell.L0*.mgikit.sample_stats`

4. `flowcell.L0*.mgikit.undetermined_barcode.complete`

This report contains the undetermined barcodes including their frequency.

5. `flowcell.L0*.mgikit.undetermined_barcode`

This report contains the top 50 frequent barcodes from the above report (4). This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

6. `flowcell.L0*.mgikit.ambiguous_barcode.complete`

This report contains the ambiguous barcodes (matched to multiple samples with the same mismatch count) including the frequency of these barcodes.

7. `flowcell.L0*.mgikit.ambiguous_barcode`

This report contains the top 50 frequent barcodes from the above report (6). This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

The first three reports must be generated for each run. It is unlikely that the fourth and fifth reports will not be generated as usually there should be some undetermined reads in the run. It is highly likely that the sixth and seventh reports will not be generated. If they are generated, it is recommended to make sure that the input sample sheet does not have issues and that the allowed mismatches are less than the minimal Hamming distance between samples.

### Execution examples

You can use the datasets at `*` to do these tests.

**Case 1:**

Running the tool and getting the help message.

```bash
./mgikit -h

```
