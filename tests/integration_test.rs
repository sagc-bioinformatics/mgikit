//use mgikit::*;
use md5;
use std::fs::File;
use std::fs;
use std::io::Read;
use std::process::Command;
use std::path::Path;
use flate2::read::MultiGzDecoder;
use walkdir::WalkDir;

fn get_hash(file_path: &String) -> Vec<u8> {
    println!("Getting hash for the file {}.", file_path);
    let mut f = File::open(file_path).unwrap();
    let mut buffer:Vec<u8> = Vec::new();
    f.read_to_end(&mut buffer).unwrap();
    buffer
}
/*
fn get_gzip_hash_old(file_path: &String) -> String {
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
*/
fn get_gzip_hash(file_path: &String) -> String {
    // changed to compare strings instead of hashs
    println!("getting string for {}", file_path);
    let mut tmp_str: String = String::new();
    let mut reader = MultiGzDecoder::new(
        File::open(Path::new(&file_path)).expect("Could not open the file")
    );
    reader.read_to_string(&mut tmp_str).unwrap();
    tmp_str
}

fn count_files_recursive(path: &String) -> u64 {
    let mut count = 0;

    for entry in WalkDir::new(path).follow_links(true) {
        match entry {
            Ok(entry) => {
                if entry.file_type().is_file() {
                    count += 1;
                }
            }
            Err(err) => eprintln!("Error while processing entry: {}", err),
        }
    }

    count
}

#[test]
fn testing_template() { 
    for ds_itr_tmp in 1..8{
        if [4, 5, 6].contains(&ds_itr_tmp){
            continue;
        }
        
        let ds_itr_in : usize = match ds_itr_tmp{
            7 =>  1,
            _ => ds_itr_tmp
         };
 
        
        let read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_1.fq.gz", ds_itr_in, ds_itr_in));
        let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_2.fq.gz", ds_itr_in, ds_itr_in));
        let out_sample_sheet_file_path : String = String::from(format!("testing_data/output/ds0{}/sample_sheet", ds_itr_tmp));
        let expected_sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.tsv", ds_itr_tmp));
        

        if ! Path::new(&format!("testing_data/output/ds0{}", ds_itr_tmp)).is_dir(){
            fs::create_dir_all(format!("testing_data/output/ds0{}", ds_itr_tmp)).unwrap();
        }
        
        let sample_sheet_file_path : String = match ds_itr_tmp{
           7 =>  String::from(format!("testing_data/input/ds0{}/sample_sheet-ds07.tsv", 1)),
           _ => String::from(format!("testing_data/input/ds0{}/sample_sheet.tsv", ds_itr_in))
        };

        //let mut digest_new;
        //let mut digest_original;
        let mut barcode_length = 0;
        let mut use_popular_template = true;
        
        if ds_itr_tmp == 2 || ds_itr_tmp == 3{
            use_popular_template = false;
        }
        if ds_itr_tmp == 3{
            barcode_length = 20;
        }
        
        
        let command = "target/debug/mgikit";
        let mut my_args: Vec<String> = vec!["template".to_string(),
                                            "-f".to_string(),
                                            read1_file_path.to_string(), 
                                            "-r".to_string(), 
                                            read2_file_path.to_string(), 
                                            "-s".to_string(), 
                                            sample_sheet_file_path.to_string(), 
                                            "--barcode-length".to_string(), 
                                            barcode_length.to_string(), 
                                            "-o".to_string(),
                                            out_sample_sheet_file_path.clone()
                                            ];
        
        if use_popular_template{
            my_args.push("--popular-template".to_string());
        }
        
        let output = Command::new(command)
            .args(my_args)
            .output() // Capture the output of the command.
            .expect("Failed to execute command");
        


        if output.status.success() {
            let output_str = String::from_utf8_lossy(&output.stdout);
            let lines: Vec<String> = output_str.split("\n").map(|it| it.to_string()).collect();
            let meta_info: Vec<String> = lines[1].split(" ").map(|it| it.to_string()).collect();
            println!("Command output:\n{} -> {}", meta_info[1], output_str);
            
        } else {
            panic!(
                "Command failed with exit code: {}\nError message: {}",
                output.status,
                String::from_utf8_lossy(&output.stderr)
            );
        }
        
        
        let digest_new = md5::compute(get_hash(&format!("{}_template.tsv", out_sample_sheet_file_path)));
        let digest_original: md5::Digest = md5::compute(get_hash(&expected_sample_sheet_file_path));
        assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
        

    }    
}

#[test]
fn testing_demultiplex() {
    for ds_itr_tmp in 1..11{
        let mut disable_illumina_format = false;
        let ds_itr_in = match ds_itr_tmp{
            6 => 1,
            9 => 8,
            7 => 1,
            4 => 3,
            5 => {disable_illumina_format = true; 1},
            _ => ds_itr_tmp
        };

        let ds_itr_ex = match ds_itr_tmp{
            6 => 1,
            _ => ds_itr_tmp
        };

        let ds_itr_fc = match ds_itr_in{
            10 => 1,
            _ => ds_itr_in
        };

        let input_folder_path = match ds_itr_tmp == 6 {
            false => String::new(),
            true => String::from(format!("testing_data/input/ds0{}/L01/", ds_itr_in))
        };
        
        let read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_1.fq.gz", ds_itr_in, ds_itr_fc));
        let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_2.fq.gz", ds_itr_in, ds_itr_fc));
        let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.tsv", ds_itr_ex));
        let lane = String::from("L01");
        let mut instrument = String::from("instrument_1"); 
        let mut run = String::from("20231212"); 
        let mut comprehensive_scan = false;
        
        if ds_itr_tmp == 5 || ds_itr_tmp == 2 {
            comprehensive_scan = true;
        }

        if ds_itr_tmp == 3 || ds_itr_tmp == 4 {
            instrument = String::from("instrument_3"); 
            run = String::from("20230727"); 
        }
    
        for allowed_mismatches in 0..5 {
            let ouput_dir = format!("testing_data/output/ds0{}/out_real-{}/", ds_itr_tmp, allowed_mismatches);
            let original_path = format!("testing_data/expected/ds0{}/ds0{}-{}/", ds_itr_ex, ds_itr_ex, allowed_mismatches);
           
            let command = "target/debug/mgikit";
            let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
                                                "-f".to_string(),
                                                read1_file_path.to_string(), 
                                                "-r".to_string(), 
                                                read2_file_path.to_string(), 
                                                "-i".to_string(), 
                                                input_folder_path.to_string(), 
                                                "-s".to_string(), 
                                                sample_sheet_file_path.to_string(), 
                                                "--lane".to_string(), 
                                                lane.to_string(), 
                                                "--run".to_string(), 
                                                run.to_string(), 
                                                "--instrument".to_string(), 
                                                instrument.to_string(), 
                                                "--writing-buffer-size".to_string(), 
                                                "131072".to_string(), 
                                                "-o".to_string(),
                                                ouput_dir.to_string(), 
                                                "-m".to_string(), 
                                                format!("{}", allowed_mismatches), 
                                                "--force".to_string()];
            
            
            if comprehensive_scan{
                my_args.push("--comprehensive-scan".to_string());
            }
            if disable_illumina_format{
                my_args.push("--disable-illumina".to_string());
            }
            if ds_itr_tmp ==  9{
                my_args.push("--template".to_string());
                my_args.push("i78:--8".to_string());
                
            }
            
            println!("{:?}", vec!["demultiplex".to_string(),
            "-f".to_string(),
            read1_file_path.to_string(), 
            "-r".to_string(), 
            read2_file_path.to_string(), 
            "-i".to_string(), 
            input_folder_path.to_string(), 
            "-s".to_string(), 
            sample_sheet_file_path.to_string(), 
            "--lane".to_string(), 
            lane.to_string(), 
            "--run".to_string(), 
            run.to_string(), 
            "--instrument".to_string(), 
            instrument.to_string(), 
            "-o".to_string(), 
            ouput_dir.to_string(), 
            "-m".to_string(), 
            format!("{}", allowed_mismatches), 
            "--force".to_string()]);


            let output = Command::new(command)
                .args(my_args)
                .output() // Capture the output of the command.
                .expect("Failed to execute command");
            
            if output.status.success() {
                let output_str = String::from_utf8_lossy(&output.stdout);
                let lines: Vec<String> = output_str.split("\n").map(|it| it.to_string()).collect();
                let meta_info: Vec<String> = lines[1].split(" ").map(|it| it.to_string()).collect();
                println!("Command output:\n{} -> {}", meta_info[1], output_str);
                
            } else {
                panic!(
                    "Command failed with exit code: {}\nError message: {}\nOutput:{}",
                    output.status,
                    String::from_utf8_lossy(&output.stderr),
                    String::from_utf8_lossy(&output.stdout)
                );
            }
            
            let paths = fs::read_dir(&original_path).unwrap();
            for path in paths {
                println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                                format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                
                if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                    
                    let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                    
                    let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                    assert_eq!(crc_new, crc_original);

                }else{
                    
                    let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                    let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                    assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
                
                }

                
        
            }

            println!("Checking count of files");
            assert_eq!(count_files_recursive(&ouput_dir),
                       count_files_recursive(&original_path));

            if [7, 8, 9, 10].contains(&ds_itr_tmp){
                break;
            }
            
        }
    }
        
}

#[test]
fn testing_reformat() {
    for mut ds_itr_tmp in 1..8{
        println!("Testing iteration: {}", ds_itr_tmp);
        let mut expected_itr: i32;
        let mut input_itr: i32;
        let mut output_itr: i32;
        
        match ds_itr_tmp{ 
            4 =>{
                expected_itr = 1;
                input_itr = 1;
                output_itr = 1;
            }, 
            5 =>{
                expected_itr = ds_itr_tmp;
                input_itr = 1;
                output_itr = ds_itr_tmp;
            },
            6 =>{
                expected_itr = ds_itr_tmp;
                input_itr = 1;
                output_itr = ds_itr_tmp;
            },
            7 =>{
                expected_itr = ds_itr_tmp;
                input_itr = 1;
                output_itr = ds_itr_tmp;
            },
            _ =>{
                expected_itr = ds_itr_tmp;
                input_itr = ds_itr_tmp;
                output_itr = ds_itr_tmp;
            }, 
            
          }

        let mut disable_illumina_format = false;
        let read1_file_path: String;
        let read2_file_path: String;
        let lane = String::from("L01");
        let mut instrument = String::from("instrument_1"); 
        let mut run = String::from("20231212"); 
        
        if ds_itr_tmp == 3{
            read1_file_path = String::from(format!("testing_data/input/extras_test/FC03_L01_sample{}_1.fq.gz", input_itr));
            read2_file_path = String::from(format!("testing_data/input/extras_test/FC03_L01_sample{}_2.fq.gz", input_itr));
            instrument = String::from("instrument_3");
            run = String::from("20230727"); 
        }else{
            read1_file_path = String::from(format!("testing_data/input/extras_test/FC01_L01_sample{}_1.fq.gz", input_itr));
            read2_file_path = String::from(format!("testing_data/input/extras_test/FC01_L01_sample{}_2.fq.gz", input_itr));
            
        }
        
        let ouput_dir     = format!("testing_data/output/extras_test/ds{}/", output_itr);
        let original_path = format!("testing_data/expected/extras_test/ds{}/", expected_itr);
           
        let command = "target/debug/mgikit";
        let mut my_args: Vec<String> = vec!["reformat".to_string(),
                                                "-f".to_string(),
                                                read1_file_path.to_string(), 
                                                "-r".to_string(), 
                                                read2_file_path.to_string(), 
                                                "--lane".to_string(), 
                                                lane.to_string(), 
                                                "-o".to_string(), 
                                                ouput_dir.to_string(),
                                                "--force".to_string(),
                                                "--sample-index".to_string(),
                                                input_itr.to_string(),];
            
        

        if disable_illumina_format{
            my_args.push("--disable-illumina".to_string());
        }
        
        if ds_itr_tmp == 2{
            my_args.push("--umi-length".to_string());
            my_args.push("8".to_string());
            
        }
        if ds_itr_tmp == 3{
            my_args.push("--umi-length".to_string());
            my_args.push("20".to_string());
        }

        if [1, 2, 3, 5].contains(&ds_itr_tmp){
            my_args.push("--run".to_string());
            my_args.push(run.to_string());
            my_args.push("--instrument".to_string());
            my_args.push(instrument.to_string());
        }else{
            my_args.push("--info-file".to_string());
            my_args.push("testing_data/input/extras_test/BioInfo.csv".to_string());
            
        }
        if ds_itr_tmp == 5{
            my_args.push("--disable-illumina".to_string());
        }
        if ds_itr_tmp == 6{
            my_args.push("--report-level".to_string());
            my_args.push("0".to_string());
        }
        if ds_itr_tmp == 7{
            my_args.push("--barcode".to_string());
            my_args.push("TCCGTTGAAA".to_string());
        }
        

        //println!("{:?}", vec![]);


        let output = Command::new(command)
            .args(my_args)
            .output() // Capture the output of the command.
            .expect("Failed to execute command");
        
        if output.status.success() {
            let output_str = String::from_utf8_lossy(&output.stdout);
            let lines: Vec<String> = output_str.split("\n").map(|it| it.to_string()).collect();
            let meta_info: Vec<String> = lines[1].split(" ").map(|it| it.to_string()).collect();
            println!("Command output:\n{} -> {}", meta_info[1], output_str);
            
        } else {
            panic!(
                "Command failed with exit code: {}\nError message: {}\nOutput:{}",
                output.status,
                String::from_utf8_lossy(&output.stderr),
                String::from_utf8_lossy(&output.stdout)
            );
        }
        
        let paths = fs::read_dir(&original_path).unwrap();
        for path in paths {
            println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                            format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
            
            if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                
                let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                
                let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                assert_eq!(crc_new, crc_original);

            }else{
                
                let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
            
            }

        }

        println!("Checking count of files");
        assert_eq!(count_files_recursive(&ouput_dir),
                    count_files_recursive(&original_path));
        
    }
}
        

