# Generate fastq files
The script `generate_fastq.py` generates large fastq files with a specific length and attaches a barcode to them for specific indices allowing few mismatches in the barcode.

## Parameters

 + **`-o or --output-file`**: Path and prefix of the of the output FASTQ files.

+ **`-n or --num-sequences`**:  Number of reads to generate (default: 1000)

+ **`-l or --sequence-length`**: Length of each sequence (default: 100)

+ **`--i7`**: i7 sequence (6 to 12 bases).

+ **`--i5`**: i5 sequence (6 to 12 bases).

+ **`--umi-len`**: UMI length should be less than 13 bases (default: 0).

+ **`--allowed-mismatches`**: Allowed mismatches, should be in [0, 1, 2]. (default: 0)

## Example

```
python3 testing_data/generate_fastq/generate_fastq.py \
    -o DS01 \
    -n 20 \
    -l 50 \
    --i7 AAAAAAAA \
    --i5 GGGGGGGG \
    --allowed-mismatches 1
```
You can create multiple samples and then merge them together using `cat` command.

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




