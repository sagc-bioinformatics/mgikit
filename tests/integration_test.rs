use mgikit::*;
use md5;
use std::fs::File;
use std::fs;
use std::io::prelude::*;


fn get_hash(file_path: &String) -> Vec<u8> {
    let mut f = File::open(file_path).unwrap();
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer).unwrap();
    buffer
}

#[test]
fn testing_ds01_mismatches() {
    let input_folder_path = String::new(); 
    let template = String::new();
    let i7_rc = false;
    let i5_rc = false;

    let disable_illumina_format = false;
    let keep_barcode = false;

    let writing_threshold = 1000;
    let read_merging_threshold = 10000;
    let comprehensive_scan = false;
    let undetermined_label = String::from("Undetermined");
    let ambiguous_label = String::from("Ambiguous");
    let force = true;
    let report_limit: usize = 50; 

    let mut read1_file_path : String = String::from("testing_data/ds01/L01/FC01_L01_R1.fastq.gz");
    let mut read2_file_path : String = String::from("testing_data/ds01/L01/FC01_L01_R2.fastq.gz");
    let sample_sheet_file_path : String = String::from("testing_data/ds01/sample_sheet.tsv");
    let lane = String::from("L01");
    let instrument = String::from("instrument_1"); 
    let run = String::from("20231212"); 
    
    let report_dir = String::new();
    
    let mut digest_new;
    let mut digest_original;
    

    for allowed_mismatches in 0..4 {
        let ouput_dir = format!("testing_data/ds01/out_real-{}/", allowed_mismatches);
        let original_path = format!("testing_data/expected_out/ds01-{}/", allowed_mismatches);
        demultiplex(
            &input_folder_path,
            &mut read2_file_path,
            &mut read1_file_path,
            &sample_sheet_file_path,
            &ouput_dir,
            &report_dir,
            allowed_mismatches,
            &template,
            i7_rc,
            i5_rc,
            &lane,
            &instrument,
            &run,
            disable_illumina_format,
            keep_barcode,
            writing_threshold,
            read_merging_threshold,
            comprehensive_scan,
            &undetermined_label,
            &ambiguous_label,
            force,
            report_limit
        );
        
        let paths = fs::read_dir(original_path).unwrap();
        for path in paths {
            println!("Checking: {}", &path.as_ref().unwrap().path().display());
            
            digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
            
            digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
            assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
    
        }
        
    }
    
    
    

    
}
