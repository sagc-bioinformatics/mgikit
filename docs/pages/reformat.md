---
title: Instructions for reformat functionality
contributors: [Ziad Al-Bkhetan]
description: User guide for MGIKIT reformat functionality including parameters details and usage examples.
toc: true
type: guides
---

## Introduction

This functionality is performed with the command `reformat`. It is to reformat reads demultiplexed by `splitBarcode` tool or raw fastq provided by MGI into illumia format and generates quality reports explained at [mgikit reports page](/mgikit/demultiplex#demultipexing-reports-section).

This command should be used for each sample separately (either paired-end or single-end). if you have multiple samples, you need to process each of them individually.

## Considerations

- If the sample barcode available in the read header (header ends with `*:N:0:********`). this barcode will be used.
- In case raw data provided and no barcode information, if the user provide sample barcode, it will be written into the header.
- If the user does not provide sample barcode or the reads don't include the sample barcode in the header, there will be no mismatches considered in the quality reports.
- If the data does not fit MGI format or splitBarcode outputs, the reformat might not be relaiable.

## Command arguments

- **`-f or --read1`**: the path to the forward reads fastq file for both paired-end and single-end input data.

- **`-r or --read2`**: the path to the reverse reads fastq file.

- **`-i or --input`**: the path to the directory that contains the input fastq files.

{% include callout.html type="note" content="Either `-i` or `-f/-r`, `-f` should be provided for a run." %}

- **`-o or --output`**: The path the output directory.

  The tool will create the directory if it does not exist
  or overwrite the content if the directory exists and the parameter `--force` is used. The tool will exit
  with an error if the directory exists, and `--force` is not used. If this parameter is not provided, the tools
  will create a directory (in the working directory) with a name based on the date and time
  of the run as follows `mgiKit_Y-m-dTHMS`. where `Y`, `m`, `d`, `H`, `M`, and `S` are the date and time format.

- **`--reports`**: The path of the output reports directory.

  By default, the tool writes the files of the run reports in the same output directory as the

  demultiplexed fastq files (`-o` or `--output` parameter). This parameter is used to write the reports in
  a different folder as specified with this parameter.

- **`--lane`**: Lane number such as `L01`.

  This parameter is used to provide the lane number when the parameter `-i` or `--input` is not

  provided. The lane number is used for QC reports and it is mandatory when Illumina format is requested for file naming.

- **`--instrument`**: The id of the sequncing machine.

  This parameter is used to provide the instrument id when the parameter `-i` or `--input` is not provided. The parameter is mandatory when Illumina format is requested for read header and file naming.

- **`--run`**: The run id. It is taken from Bioinf.csv as the date and time of starting the run.

  This parameter is used to provide the run id when the parameter `-i` or `--input` is not provided. The parameter is mandatory when Illumina format is requested for read header and file naming.

- **`--writing-buffer-size`**: The default value is `67108864`. The size of the buffer for each sample to be filled with data then written once to the disk. Smaller buffers will need less memory but makes the tool slower. Largeer buffers need more memory.

- **`--compression-level`**: The level of compression (between 0 and 12). 0 is fast but no compression, 12 is slow but high compression. [default: 1]

- **`--force`**: this flag is to force the run and overwrite the existing output directory if exists.

- **`--info-file`**: The name of the info file that contains the run information. Only needed when using the `--input` parameter. [default: BioInfo.csv]

- **`--disable-illumina`**: reads will be left as is and only quality reports will be generated.

- **`--umi-length`**: The length of UMI expected at the end of the read (r1 for single-end, or r2 for paired-end) [Default: 0].

- **`--report-level`**: The level of reporting. 0 no reports will be generated, 1 data quality and demultiplexing reports. 2: all reports (reports on data quality, demultiplexing, undetermined and ambiguous barcodes).[default: 2]

- **`--sample-index`**: The index of the sample in the sample sheet. It is required for file naming. [default: 1]

- **`--barcode`**: The barcode of the specific sample to calculate the mismatches for the reports. If not provided, no mismatches will be calculated.

- **`--validate`**: when enabled, the tool will validate the content of the input fastq files.

## Usage Examples

**1. Demultiplexing a run with dual indexes (i7 and i5)**

```bash
target/release/mgikit  reformat \
    -f testing_data/input/extras_test/FC01_L01_sample1_1.fq.gz \
    -r testing_data/input/extras_test/FC01_L01_sample1_2.fq.gz \
    --lane L01 -o output \
    --sample-index 1 \
    --info-file testing_data/input/extras_test/BioInfo.csv
```
