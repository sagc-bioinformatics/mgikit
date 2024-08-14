---
title: Instructions for template detection
contributors: [Ziad Al-Bkhetan]
description: User guide for MGIKIT template detection functionality including parameters details and usage examples.
toc: true
type: guides
---

This command is used to detect the location and form of the indexes within the read barcode. It simply goes through a small number of the reads and investigates the number of matches with the indexes in the sample sheet within each possible location in the read barcode and considers the indexes as is and their reverse complementary.

It reports matches for all possible combinations and uses the read template that had the maximum number of matches. This process happens for each sample individually and therefore, the best matching template for each sample will be reported.

Using this comprehensive scan, the tool can detect the templates for mixed libraries.

## Parameters

**Fastq input file**

- **`-f or --read1`**: the path to the forward reads fastq file for both paired-end and single-end input data.

- **`-r or --read2`**: the path to the reverse reads fastq file.

- **`-s or --sample-sheet`**: the path to the sample sheet file.

  This is the same format as above, but only sample_id and i7 are required. i5 is required for dual indexes data.

- **`-o or --output`**: The path and prefix of output files. The tools will create two files at the same path with the same prefix and end with `_template.tsv` and `_details.tsv`.

- **`--testing-reads`**: The number of reads to be investigated to check and detect the templates. The default is 5,000 reads. A Larger number increases the performance time.

- **`--barcode-length`**: The length of the read barcode at the end of the read2 in paired-end or read1 in single end to be investigated. By default, the barcode length is set to be the length difference between read2 and read1.

- **`--no-umi`**: If the barcode contains extra base pairs other than the indexes, the tool considers the longest as an umi. If this parameter is enabled, the tool will ignore all extra base pairs in the barcode and trim them from the read.

- **`--popular-template`**: by default, the tool reports the template that matches the maximum number of reads to each corresponding sample. If this option is enabled, the tool will use the most frequent template across all samples as the final template for all samples.

- **`--max-umi-length`**: if barcode length is not provided, the tool will set the barcode length to the length difference between read2 and read1. If the barcode length is greater than the sum of indexes lengths and this parameter, the tool will stop. The default is 10 bp. You can disable this parameter by either providing a large number or providing the barcode length (`--barcode-length`) parameter manually.
