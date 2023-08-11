use mgikit::*;
use md5;
use std::fs::File;
use std::fs;
use std::io::Read;
use std::process::Command;
use std::path::Path;

fn get_hash(file_path: &String) -> Vec<u8> {
    let mut f = File::open(file_path).unwrap();
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer).unwrap();
    buffer
}

fn get_gzip_hash(file_path: &String) -> String {
    println!("getting hash for {}", file_path);
    let command = "gzip";
    let output = Command::new(command)
                    .arg("-lv")
                    .arg(file_path)
                    .output() // Capture the output of the command.
                    .expect("Failed to execute command");

    if output.status.success() {
        let output_str = String::from_utf8_lossy(&output.stdout);
        let lines: Vec<String> = output_str.split("\n").map(|it| it.to_string()).collect();
        let meta_info: Vec<String> = lines[1].split(" ").map(|it| it.to_string()).collect();
        //println!("Command output:\n{} -> {}", meta_info[1], output_str);
        meta_info[1].trim().to_string()
    } else {
        panic!(
            "Command failed with exit code: {}\nError message: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }
}

#[test]
fn testing_template() { 
    for ds_itr in 1..4{
        let read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_R1.fastq.gz", ds_itr, ds_itr));
        let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_R2.fastq.gz", ds_itr, ds_itr));
        let sample_sheet_file_path : String = String::from(format!("testing_data/input/ds0{}/sample_sheet.tsv", ds_itr));
        let out_sample_sheet_file_path : String = String::from(format!("testing_data/output/ds0{}/sample_sheet", ds_itr));
        let expected_sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.tsv", ds_itr));
        
        if ! Path::new(&format!("testing_data/output/ds0{}", ds_itr)).is_dir(){
            fs::create_dir_all(format!("testing_data/output/ds0{}", ds_itr)).unwrap();
        }
        

        //let mut digest_new;
        //let mut digest_original;
        let mut barcode_length = 0;
        let mut use_popular_template = true;
        
        if ds_itr == 2 || ds_itr == 3{
            use_popular_template = false;
        }
        if ds_itr == 3{
            barcode_length = 20;
        }
        detect_template(&read1_file_path, 
            &read2_file_path,
            &sample_sheet_file_path,
            &out_sample_sheet_file_path,
            1000,
            barcode_length,
            true,
            use_popular_template,
            10
        );
        
        let digest_new = md5::compute(get_hash(&format!("{}_template.tsv", out_sample_sheet_file_path)));
        let digest_original: md5::Digest = md5::compute(get_hash(&expected_sample_sheet_file_path));
        assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
        

    }    
}

#[test]
fn testing_demultiplex() {
    for ds_itr_tmp in 1..6{
        
        let mut disable_illumina_format = false;
        let ds_itr = match ds_itr_tmp{
            4 => 3,
            5 => {disable_illumina_format = true; 1},
            _ => ds_itr_tmp
        };

        let input_folder_path = String::new(); 
        let template = String::new();
        let i7_rc = false;
        let i5_rc = false;
    
        let read1_file_name_suf: String =  String::from("_read_1.fq.gz");
        let read2_file_name_suf: String =  String::from("_read_2.fq.gz");
        let info_file: String =  String::from("BioInfo.csv");
        
    
        let keep_barcode = false;
    
        let writing_threshold = 1000;
        let read_merging_threshold = 10000;
        let mut comprehensive_scan = false;
        let undetermined_label = String::from("Undetermined");
        let ambiguous_label = String::from("Ambiguous");
        let force = true;
        let report_limit: usize = 50; 
    
        let mut read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_R1.fastq.gz", ds_itr, ds_itr));
        let mut read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_R2.fastq.gz", ds_itr, ds_itr));
        let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.tsv", ds_itr_tmp));
        let lane = String::from("L01");
        let mut instrument = String::from("instrument_1"); 
        let mut run = String::from("20231212"); 
        
        let report_dir = String::new();
        
        let mut digest_new;
        let mut digest_original;
        
        if ds_itr_tmp == 5 || ds_itr_tmp == 2 {
            comprehensive_scan = true;
        }

        if ds_itr_tmp == 3 || ds_itr_tmp == 4 {
            instrument = String::from("instrument_3"); 
            run = String::from("20230727"); 
        }
    
        for allowed_mismatches in 0..5 {
            let ouput_dir = format!("testing_data/output/ds0{}/out_real-{}/", ds_itr_tmp, allowed_mismatches);
            let original_path = format!("testing_data/expected/ds0{}/ds0{}-{}/", ds_itr_tmp, ds_itr_tmp, allowed_mismatches);
            demultiplex(
                &input_folder_path,
                &mut read1_file_path,
                &mut read2_file_path,
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
                report_limit,
                &read1_file_name_suf,
                &read2_file_name_suf,
                &info_file
            );
            
            let paths = fs::read_dir(original_path).unwrap();
            for path in paths {
                println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                                format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap().clone()));
                
                if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                    
                    let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap().clone()));
                    
                    let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                    assert_eq!(crc_new, crc_original);

                }else{
                    digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                    digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                    assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
                
                }

                
        
            }
            

            
        }
    }
        
}


