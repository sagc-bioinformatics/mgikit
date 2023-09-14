---

![SAGC-Bioinformatics](docs/assets/SAGC-logo-hover.png)

---
# MGIKIT 
mgikit is a collection of tools used to demultiplex fastq files and generate reports including.

The toolkit includes the following functionalities:

1. **demultiplex**: This command is used to demultiplex fastq files and assign the sequencing reads to their
associated samples.
2. **template**: This command is used to detect the location and form of the indexes within the read barcode if unknown.
3. **reports**: This command is used to merge Lanes' reports (demultiplexing and quality reports) and generate the same reports for the whole run and for grouped samples (multiple projects in the same run).

Details about these functionalities are described at [mgikit usage documentation](docs/usage_documentation.md). 

# Installation

The binary file can be used directly, or the tool can be built from the source code as follows:

```bash
clone 
cd source_code_path
cargo build --release
```

# Note

This is my first project using the Rust programming language. Very likely that the code can be improved and optimised. I will be improving it with time.   

# Commerical Use

Please contact us if you want to use the software for commercial purposes.

# Contact information

For any help or inquiries, please contact: ziad.albkhetan@gmail.com
