---
title: Instructions to generate MultQC report from MGIKIT reports.
contributors: [Ziad Al-Bkhetan]
description: User guide for generating html report summarising the demultiplexing results using mgikit output reports.
toc: true
type: guides
---

## mgikit Reports

The demultiplex command generates multiple reports with file names that start with the flowcell and lane being demultiplexed.
a MultiQC hitm report can be generated from these reports using [mgikit-multiqc](https://github.com/sagc-bioinformatics/mgikit-multiqc) plugin as described at the plugin [repository](https://github.com/sagc-bioinformatics/mgikit-multiqc).

1. `flowcell.L0*.mgikit.info`

This report contains the number of reads per sample respectively to each possible mismatch.
For example:

| **sample** | **0-mismatches** | **1-mismatches** |
| :--------: | :--------------: | :--------------: |
|    S01     |       3404       |       5655       |

This means that there was only one mismatch allowed during this execution and the sample S01 has 3404 reads with indexes matching perfectly and 5655 reads with indexes that differ by 1 base compared to the indexes provided in the sample sheet.

This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

2. `flowcell.L0*.mgikit.general`

This file contains summary information related to the cluster count and quality scores, summarised for each sample as well as at the whole lane scale. This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

**Report content:**

- **Lane statistics columns**

  1. `Run ID-Lane`: Run ID and lane number.
  2. `Mb Total Yield`: total number of bases in a million.
  3. `M Total Clusters`: total number of reads in million.
  4. `% bases ≥ Q30`: percentage of bases with a quality score greater than 30 of all bases.
  5. `Mean Quality`: The average quality score for the bases.
  6. `% Perfect Index`: The percentage of reads with perfectly matching indices of all reads.

- **Sample general info**

  1. `Sample ID`: sample ID taken from the sample sheet.
  2. `M Clusters`: total number of reads in million.
  3. `Mb Yield ≥ Q30`: total number of bases with a quality score greater than 30 in million.
  4. `% R1 Yield ≥ Q30`: percentage of bases with a quality score greater than 30 of all bases calculated only for forward reads.
  5. `% R2 Yield ≥ Q30`: percentage of bases with a quality score greater than 30 of all bases calculated only for reverse reads.
  6. `% R3 Yield ≥ Q30`: percentage of bases with a quality score greater than 30 of all bases calculated only for indices.
  7. `% Perfect Index`: The percentage of reads with perfectly matching indices of all reads.

3. `flowcell.L0*.mgikit.sample_stats`

This file contains the informaiton in the above mentioned reports but in simple format. This is used to merge the reports from multiple lanes into one report for the whole run.

4. `flowcell.L0*.mgikit.undetermined_barcode.complete`

This report contains the undetermined barcodes including their frequency.

5. `flowcell.L0*.mgikit.undetermined_barcode`

This report contains the top 50 frequent barcodes from the above report (4). This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

6. `flowcell.L0*.mgikit.ambiguous_barcode.complete`

This report contains the ambiguous barcodes (matched to multiple samples with the same mismatch count) including the frequency of these barcodes.

7. `flowcell.L0*.mgikit.ambiguous_barcode`

This report contains the top 50 frequent barcodes from the above report (6). This file is used for the mgikit plugin to visualise quality control reports through MultiQC.

The first three reports must be generated for each run. It is unlikely that the fourth and fifth reports will not be generated as usually there should be some undetermined reads in the run. It is highly likely that the sixth and seventh reports will not be generated. If they are generated, it is recommended to make sure that the input sample sheet does not have issues and that the allowed mismatches are less than the minimal Hamming distance between samples.

## Example: Generate MultiQC report from mgikit reports

In order to generate a [multiqc](https://multiqc.info/) report from mgikit reports, multiqc needs to be installed.

Here is an example of how to generate the report:

```bash
pip install MultiQC
git clone https://github.com/sagc-bioinformatics/mgikit-multiqc

cd mgikit-multiqc
python setup.py install

# To test the installation, run the following on the testing dataset available at the plugin repository
multiqc mgikit-examples/test/

```
