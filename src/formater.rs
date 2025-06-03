use getset::{ Getters, CopyGetters };
use log::{ info, warn };

#[derive(Getters, CopyGetters, Clone, Default)]
pub struct ReformatedSample {
    #[getset(get_copy = "pub")]
    sample_index: usize,
    #[getset(get_copy = "pub")]
    umi_length: usize,
    #[getset(get = "pub")]
    barcode: String,
    #[getset(get = "pub")]
    sample_label: String,
}

impl ReformatedSample {
    pub fn new(
        sample_label: String,
        sample_index: usize,
        umi_length: usize,
        barcode: String
    ) -> Self {
        if sample_index < 1 {
            panic!("Sample index (val: {}) needs to be greater than 0!", sample_index);
        }
        info!("Sample's id, extracted from the file name is: {}", sample_label);
        info!("Sample's index to be used in file naming is: {}", sample_index);
        info!("Sample's barcode: {}", barcode);
        if barcode.len() == 0 {
            warn!(
                "No mismatches will be considered as sample's barcode is not provided to compare against!"
            );
        }
        info!("Umi length to be extracted from the tail of R2: {}", umi_length);
        Self {
            sample_label,
            sample_index,
            umi_length,
            barcode,
        }
    }
}

pub fn parse_sb_file_name(file_name: &String) -> (String, String, String, String) {
    //V350170513_L01_23-03891_1.fq.gz
    let parts: Vec<&str> = file_name.split("_").collect();
    if parts.len() < 3 {
        warn!(
            "Input file name is expected to have the following format 'Flowcell_Lane_SampleLabel_1/2.fq/fastq.gz'"
        );
        warn!("Ignore this warning if you provided Sample label and lane parameters.");
        (String::new(), String::new(), String::new(), String::new())
    } else {
        (
            parts[2..parts.len() - 1].join("_"),
            parts[0].to_string(),
            parts[1].to_string(),
            parts[parts.len() - 1].to_string(),
        )
    }
}
