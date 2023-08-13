---
title: Instructions for report functionality
contributors: [Ziad Al-Bkhetan]
description: User guide for MGIKIT report functionality including parameters details and usage examples.
toc: true
type: guides
---

## Introduction

mgikit demultiplexer outputs four quality information and demultiplexing reports that can be used for reporting and can be read by MultiQC tool and mgikit plugin.
These reports are explained at [mgikit reports page](/mgikit/demultiplex#demultipexing-reports-section).

if the run has multiple lanes, there will be lane-specific reports. The reports can be for the whole run and each project within this run. The lanes generated from each lane can be merged to generate quality and demultiplexing reports for all samples by merging the information from the lanes' reports. This `report` command does this merge for you.

## Command arguments

+ **`--qc-reports`**: The paths to the QC reports are separated by a space ` `. This argument takes multiple values and is mandatory. The tool expects here the reports generated for each lane in the run and you also can combine the reports generated from multiple runs for the same samples.
 

+ **`-o or --output`**: The path and prefix of output files. The tools will create two files at the same path with the same prefix and end with `.info` and `.general`.

## Usage Examples


