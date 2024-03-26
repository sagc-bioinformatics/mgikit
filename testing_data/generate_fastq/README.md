# Generate DS03 dataset
Run the following command
```
python3 generate_ds03.py
```
This will create a dataset with paired-end files `DS03_R1.fastq.gz` and `DS03_R2.fastq.gz` in the working directory.
You can change the path and dataset name by editing the script (`output_dir`).

# Demultipex with mgikit
Run the following command
```
mgikit demultiplex -f DS03_R1.fastq.gz -r DS03_R2.fastq.gz -s generate_fastq/ds03_samplesheet_mk.tsv -o ./mgikit_output --lane L01 --disable-illumina
```

# Demultipex with splitBarcode
Run the following command
```
splitBarcode -m 7 -t 1 -1 DS03_R1.fastq.gz -2  DS03_R2.fastq.gz -o -B generate_fastq/ds03_samplesheet_sb.tsv -b 200 8 1 -b 216 8 1 -o ./sb_output
```




