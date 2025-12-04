//use mgikit::*;
use md5;
use std::collections::HashMap;
use std::fs::File;
use std::fs;
use std::io::Read;
use std::process::Command;
use std::path::{Path, PathBuf};
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
        
        print!("Command: {:?}", &my_args.join(" "));
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
    for threads_cnt in [1, 2, 3, 4, 5, 8, 10]{
        for ds_itr_tmp in 1..15{
                let mut disable_illumina_format = false;
                let ds_itr_in = match ds_itr_tmp{
                    6 => 1,
                    9 => 8,
                    7 => 1,
                    4 => 3,
                    5 => {disable_illumina_format = true; 1},
                    11 => 1,
                    12 => 2,
                    14 => 2,
                    _ => ds_itr_tmp
                };

                let ds_itr_ex = match ds_itr_tmp{
                    6 => 1,
                    14 => 12,
                    _ => ds_itr_tmp
                };

                let ds_itr_fc = match ds_itr_in{
                    10 => 1,
                    11 => 1,
                    12 => 1,
                    13 => 2,
                    14 => 2,
                    _ => ds_itr_in
                };

                let input_folder_path = match ds_itr_tmp == 6 {
                    false => String::new(),
                    true => String::from(format!("testing_data/input/ds0{}/L01/", ds_itr_in))
                };
                
                let read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_1.fq.gz", ds_itr_in, ds_itr_fc));
                let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_2.fq.gz", ds_itr_in, ds_itr_fc));
                let ext = if ds_itr_tmp == 1 {"csv"} else {"tsv"};
                let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.{ext}", ds_itr_ex));
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
                    if PathBuf::from(&ouput_dir).exists() {
                        fs::remove_dir_all(&ouput_dir).unwrap();
                    }
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
                                                        "--force".to_string(),
                                                        "-t".to_string(),
                                                        threads_cnt.to_string()];
                                
                    if comprehensive_scan{
                        my_args.push("--comprehensive-scan".to_string());
                    }
                    my_args.push("--validate".to_string());
                    if disable_illumina_format{
                        my_args.push("--disable-illumina".to_string());
                    }
                    if ds_itr_tmp ==  9 {
                        my_args.push("--template".to_string());
                        my_args.push("i78:--8".to_string());
                        
                    }

                    if ds_itr_tmp <  11{
                        my_args.push("--all-index-error".to_string());                
                    }
                    println!("{:?}", &my_args);


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
                    if ds_itr_tmp > 10 && allowed_mismatches == 2{
                        break;
                    }
                    
                }
            }
            
    }
}

#[test]
fn testing_demultiplex_not_mgi_input() {
    //let mut disable_illumina_format = false;
    let read1_file_path : String = String::from("testing_data/input/ds015/L01/FC02_L01_read_1.fq.gz");
    let read2_file_path : String = String::from("testing_data/input/ds015/L01/FC02_L01_read_2.fq.gz");
    let sample_sheet_file_path : String = String::from("testing_data/expected/ds015/sample_sheet_expected.tsv");
    let lane = String::from("L01");
    for allowed_mismatches in 0..3 {
        let ouput_dir = format!("testing_data/output/ds015/out_real-{}/", allowed_mismatches);
        let original_path = format!("testing_data/expected/ds015/ds015-{}/", allowed_mismatches);
        if PathBuf::from(&ouput_dir).exists() {
            fs::remove_dir_all(&ouput_dir).unwrap();
        } 

        let command = "target/debug/mgikit";
        let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
                                                "-f".to_string(),
                                                read1_file_path.to_string(), 
                                                "-r".to_string(), 
                                                read2_file_path.to_string(), 
                                                "-s".to_string(), 
                                                sample_sheet_file_path.to_string(), 
                                                "--lane".to_string(), 
                                                lane.to_string(), 
                                                "--writing-buffer-size".to_string(), 
                                                "131072".to_string(), 
                                                "-o".to_string(),
                                                ouput_dir.to_string(), 
                                                "-m".to_string(), 
                                                format!("{}", allowed_mismatches), 
                                                "--force".to_string()];
                        
        my_args.push("--disable-illumina".to_string());
        my_args.push("--not-mgi".to_string());                
        my_args.push("--validate".to_string());
        println!("{:?}", &my_args);

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
        
        let mut output_paths = HashMap::new();
        for path in fs::read_dir(&ouput_dir).unwrap() {
            if path.as_ref().unwrap().file_name().to_str().unwrap().contains("L01.mgikit.general"){
                output_paths.insert("L01.mgikit.general",path.as_ref().unwrap().file_name().into_string().unwrap());
            }else if path.as_ref().unwrap().file_name().to_str().unwrap().contains("L01.mgikit.info"){
                output_paths.insert("L01.mgikit.info", path.as_ref().unwrap().file_name().into_string().unwrap());
            }else if path.as_ref().unwrap().file_name().to_str().unwrap().contains("L01.mgikit.sample_stats"){
                output_paths.insert("L01.mgikit.sample_stats", path.as_ref().unwrap().file_name().into_string().unwrap());
            }else if path.as_ref().unwrap().file_name().to_str().unwrap().contains("L01.mgikit.undetermined_barcode.complete"){
                output_paths.insert("L01.mgikit.undetermined_barcode.complete", path.as_ref().unwrap().file_name().into_string().unwrap());
            }else if path.as_ref().unwrap().file_name().to_str().unwrap().contains("L01.mgikit.undetermined_barcode"){
                output_paths.insert("L01.mgikit.undetermined_barcode", path.as_ref().unwrap().file_name().into_string().unwrap());
            }
        }

        println!("{:?}", output_paths);
        let paths = fs::read_dir(&original_path).unwrap();
        for path in paths {
            
            
            if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                println!("Checking: {} and {}", path.as_ref().unwrap().path().display(), 
                          format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));

                let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                
                let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                assert_eq!(crc_new, crc_original);

            }else{
                println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                                format!("{}{}", ouput_dir, output_paths.get(&path.as_ref().unwrap().file_name().to_str().unwrap()).unwrap()));    
                    
                let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, output_paths.get(&path.as_ref().unwrap().file_name().to_str().unwrap()).unwrap())));
                let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
            
            }

            
    
        }

        println!("Checking count of files");
        assert_eq!(count_files_recursive(&ouput_dir),
                    count_files_recursive(&original_path) + 1);        
    }
            
}

#[test]
fn testing_reformat() {
    for ds_itr_tmp in 1..8{
        println!("Testing iteration: {}", ds_itr_tmp);
        let expected_itr: i32;
        let input_itr: i32;
        let output_itr: i32;
        
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

        let disable_illumina_format = false;
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

        print!("Command: {:?}", &my_args.join(" "));
        
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
        
#[test]
fn testing_reformat_raw() {
    for ds_itr_tmp in 1..13{
        if ds_itr_tmp == 8{
            continue;
        }
        println!("Testing iteration: {}", ds_itr_tmp);
        let expected_itr: i32;
        let input_itr: i32;
        let output_itr: i32;
        
        match ds_itr_tmp{ 
            4 =>{
                expected_itr = 1;
                input_itr = 1;
                output_itr = ds_itr_tmp;
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
            12 =>{
                expected_itr = ds_itr_tmp;
                input_itr = 1;
                output_itr = 12;
            },
            v if v > 8 && v < 12 => {
                expected_itr = ds_itr_tmp;
                input_itr = v - 8;
                output_itr = ds_itr_tmp;
            },
            _ =>{
                expected_itr = ds_itr_tmp;
                input_itr = ds_itr_tmp;
                output_itr = ds_itr_tmp;
            }, 
            
          }


        let disable_illumina_format = false;
        let read1_file_path: String;
        let read2_file_path: String;
        let lane = String::from("L01");
        let mut instrument = String::from("instrument_1"); 
        let mut run = String::from("20231212"); 
        
        if ds_itr_tmp == 3 || ds_itr_tmp == 11{
            read1_file_path = String::from(format!("testing_data/input/extras_test/raw/FC03_L01_sample{}_1.fq.gz", input_itr));
            read2_file_path = String::from(format!("testing_data/input/extras_test/raw/FC03_L01_sample{}_2.fq.gz", input_itr));
            instrument = String::from("instrument_3");
            run = String::from("20230727"); 
        }else{
            read1_file_path = String::from(format!("testing_data/input/extras_test/raw/FC01_L01_sample{}_1.fq.gz", input_itr));
            read2_file_path = String::from(format!("testing_data/input/extras_test/raw/FC01_L01_sample{}_2.fq.gz", input_itr));
            
        }
        
        let ouput_dir     = format!("testing_data/output/extras_test/raw-ds{}/", output_itr);
        let original_path = format!("testing_data/expected/extras_test/ds{}{}/", expected_itr, if ds_itr_tmp == 7 {"-raw"} else {""});
           
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
        
        if ds_itr_tmp == 2 || ds_itr_tmp == 10{
            my_args.push("--umi-length".to_string());
            my_args.push("8".to_string());
            
        }
        if ds_itr_tmp == 3 || ds_itr_tmp == 11{
            my_args.push("--umi-length".to_string());
            my_args.push("20".to_string());
        }

        if [1, 2, 3, 5, 11, 12].contains(&ds_itr_tmp){
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
        
        if [1, 2, 4, 6].contains(&ds_itr_tmp){
            my_args.push("--barcode".to_string());
            my_args.push("TCCGTTGAAT".to_string());
        }else if ds_itr_tmp == 3{
            my_args.push("--barcode".to_string());
            my_args.push("AAAAAAAA".to_string());
        }
        if ds_itr_tmp == 12{
            my_args.push("--sample-label".to_string());
            my_args.push("test-sample".to_string());
        }
        //println!("{:?}", vec![]);

        print!("Command: {:?}", &my_args.join(" "));
        
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
            
            let out_file = if ds_itr_tmp == 12{
                format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap().replace("sample1", "test-sample"))
            }else{
                format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())
            };
            println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                            out_file);
            
            if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                
                let crc_new = get_gzip_hash(&out_file);
                
                let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                assert_eq!(crc_new, crc_original);

            }else{
                
                let digest_new = md5::compute(get_hash(&out_file));
                let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
            
            }

        }

        println!("Checking count of files");
        assert_eq!(count_files_recursive(&ouput_dir),
                    count_files_recursive(&original_path));
        
    }
}

#[test]
fn testing_demultiplex_threads() {
    for ds_itr_tmp in 1..15{
        let mut disable_illumina_format = false;
        let ds_itr_in = match ds_itr_tmp{
            6 => 1,
            9 => 8,
            7 => 1,
            4 => 3,
            5 => {disable_illumina_format = true; 1},
            11 => 1,
            12 => 2,
            14 => 2,
            _ => ds_itr_tmp
        };

        let ds_itr_ex = match ds_itr_tmp{
            6 => 1,
            14 => 12,
            _ => ds_itr_tmp
        };

        let ds_itr_fc = match ds_itr_in{
            10 => 1,
            11 => 1,
            12 => 1,
            13 => 2,
            14 => 2,
            _ => ds_itr_in
        };

        let ds_threads = match ds_itr_tmp{
            1 => 1,
            2 => -1,
            6 => 0,
            7 => 2,
            8 => 3,
            9 => 2,
            _ => ds_itr_tmp
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
            
            let ouput_dir = format!("testing_data/output_thread/ds0{}/out_real-{}/", ds_itr_tmp, allowed_mismatches);
            let original_path = format!("testing_data/expected/ds0{}/ds0{}-{}/", ds_itr_ex, ds_itr_ex, allowed_mismatches);
            if PathBuf::from(&ouput_dir).exists() {
                fs::remove_dir_all(&ouput_dir).unwrap();
            }
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
            if ds_itr_tmp ==  9 {
                my_args.push("--template".to_string());
                my_args.push("i78:--8".to_string());
                
            }

            if ds_threads !=  -1 {
                my_args.push("--threads".to_string());
                my_args.push(ds_threads.to_string());
                
            }
            my_args.push("--validate".to_string());
            if ds_itr_tmp <  11{
                my_args.push("--all-index-error".to_string());                
            }
            println!("{:?}", my_args);

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
            if ds_itr_tmp > 10 && allowed_mismatches == 2{
                break;
            }
            
        }
    }
        
}

#[test]
fn testing_demultiplex_large() {
    for se in 0..2{
        for thread_cnt in [1, 2, 3, 4, 5, 8, 10]{
                let read1_file_path : String = String::from(format!("testing_data/input/large_ds/ZFC01_L01_read_1.fq.gz"));
                let read2_file_path : String = String::from(format!("testing_data/input/large_ds/ZFC01_L01_read_2.fq.gz"));
                let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds01/sample_sheet_expected.tsv"));
                let lane = String::from("L01");
                let instrument = String::from("instrument_1"); 
                let run = String::from("20231212");
                
                    
                let ouput_dir = String::from("testing_data/output_large/");
                let original_path = String::from("testing_data/expected/large_ds/");
                if PathBuf::from(&ouput_dir).exists() {
                    fs::remove_dir_all(&ouput_dir).unwrap();
                }
                let command = "target/debug/mgikit";
                let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
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
                                                    format!("{}", 0), 
                                                    "--force".to_string()];
                            
                my_args.push("--threads".to_string());
                my_args.push(thread_cnt.to_string());
                my_args.push("--validate".to_string());

                my_args.push("-f".to_string());
                if se == 0{
                    my_args.push(read1_file_path.to_string());
                    my_args.push("-r".to_string());
                    my_args.push(read2_file_path.to_string());     
                }else{
                    my_args.push(read2_file_path.to_string());
                }
                
                my_args.push("--all-index-error".to_string());                
                
                println!("{:?}", my_args);

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
                if se == 0{
                    let paths = fs::read_dir(&original_path).unwrap();
                    for path in paths {
                        println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                                        format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                        if format!("{}", path.as_ref().unwrap().path().display()).contains("se-FC01.L01.mgikit"){
                            continue;
                        }

                        if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                            if thread_cnt > 1{
                                continue;
                            }
                            let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
                            
                            let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                            assert_eq!(crc_new, crc_original);

                        }else{   
                            let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                            let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                            assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
                        }
                    }       
                }
                
            }
    }
}

#[test]
fn testing_demultiplex_large_se() {
    for thread_cnt in 1..4{
        let read2_file_path : String = String::from(format!("testing_data/input/large_ds/ZFC01_L01_read_2.fq.gz"));
        let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds01/sample_sheet_expected.tsv"));
        let lane = String::from("L01");
        let instrument = String::from("instrument_1"); 
        let run = String::from("20231212"); 
        
            
        let ouput_dir = String::from("testing_data/output_large-se/");
        let original_path = String::from("testing_data/expected/large_ds/");
        if PathBuf::from(&ouput_dir).exists() {
            fs::remove_dir_all(&ouput_dir).unwrap();
        }
        let command = "target/debug/mgikit";
        let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
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
                                            format!("{}", 0), 
                                            "--force".to_string()];
                    
        my_args.push("--threads".to_string());
        my_args.push(thread_cnt.to_string());
        my_args.push("--validate".to_string());

        my_args.push("-f".to_string());
        my_args.push(read2_file_path.to_string());
    
        

        my_args.push("--all-index-error".to_string());                
        
        println!("{:?}", my_args);

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
            if format!("{}", path.as_ref().unwrap().path().display()).contains("se-FC01.L01.mgikit"){
                continue;
            }
            println!("Checking: {} and {}", path.as_ref().unwrap().path().display(),
                                            format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap()));
            
            if path.as_ref().unwrap().path().display().to_string().contains("_R1_001"){
                continue; 
            }

            if path.as_ref().unwrap().path().display().to_string().contains("sample_stats"){
                let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                let digest_original = md5::compute(get_hash(&String::from("testing_data/expected/large_ds/se-FC01.L01.mgikit.sample_stats")));
                assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
                continue; 
            }
            
            if format!("{}", &path.as_ref().unwrap().path().display()).ends_with(".gz"){
                if thread_cnt > 1{
                    continue;
                }
                let tmp = &path.as_ref().unwrap().file_name().into_string().unwrap();
                let crc_new = get_gzip_hash(&format!("{}{}", ouput_dir, &tmp.replace("R2", "R1")));
                
                let crc_original = get_gzip_hash(&format!("{}", &path.unwrap().path().display()));
                assert_eq!(crc_new, crc_original);

            }else{   
                let digest_new = md5::compute(get_hash(&format!("{}{}", ouput_dir, &path.as_ref().unwrap().file_name().to_str().unwrap())));
                let digest_original = md5::compute(get_hash(&format!("{}", &path.unwrap().path().display())));
                assert_eq!(format!("{:x}", digest_new), format!("{:x}", digest_original));
            }
        }

        
    }    
}

#[test]
fn testing_single_end() {
    for thread in 1..6{    
        for ds_itr_tmp in 1..15{
            let mut disable_illumina_format = false;
            let ds_itr_in = match ds_itr_tmp{
                6 => 1,
                9 => 8,
                7 => 1,
                4 => 3,
                5 => {disable_illumina_format = true; 1},
                11 => 1,
                12 => 2,
                14 => 2,
                _ => ds_itr_tmp
            };

            let ds_itr_ex = match ds_itr_tmp{
                6 => 1,
                14 => 12,
                _ => ds_itr_tmp
            };

            let ds_itr_fc = match ds_itr_in{
                10 => 1,
                11 => 1,
                12 => 1,
                13 => 2,
                14 => 2,
                _ => ds_itr_in
            };

            let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_2.fq.gz", ds_itr_in, ds_itr_fc));
            let ext = if ds_itr_tmp == 1 {"csv"} else {"tsv"};
            let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.{ext}", ds_itr_ex));
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
        
            for allowed_mismatches in 0..2 {
                if allowed_mismatches == 1 &&  ds_itr_tmp != 2{
                    continue;
                }
                let ouput_dir = format!("testing_data/output-se/ds0{}/out_real-{}/", ds_itr_tmp, allowed_mismatches);
                if PathBuf::from(&ouput_dir).exists() {
                    fs::remove_dir_all(&ouput_dir).unwrap();
                }
                let original_path = format!("testing_data/expected/se-test/ds0{}/ds0{}-{}/", ds_itr_ex, ds_itr_ex, allowed_mismatches);
            
                let command = "target/debug/mgikit";
                let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
                                                    "-f".to_string(),
                                                    read2_file_path.to_string(), 
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
                                                    "--force".to_string(),
                                                    "-t".to_string(),
                                                    thread.to_string()];
                            
                if comprehensive_scan{
                    my_args.push("--comprehensive-scan".to_string());
                }
                if disable_illumina_format{
                    my_args.push("--disable-illumina".to_string());
                }
                if ds_itr_tmp ==  9 {
                    my_args.push("--template".to_string());
                    my_args.push("i78:--8".to_string());
                    
                }
                my_args.push("--validate".to_string());

                if ds_itr_tmp <  11{
                    my_args.push("--all-index-error".to_string());                
                }
                
                println!("{:?}", &my_args);


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

                if ds_itr_tmp > 1{
                    continue;
                }

                println!("Checking count of files");
                assert_eq!(count_files_recursive(&ouput_dir),
                        count_files_recursive(&original_path));

                if [7, 8, 9, 10].contains(&ds_itr_tmp){
                    break;
                }
                if ds_itr_tmp > 10 && allowed_mismatches == 2{
                    break;
                }
                
            }
        }
    }
}

#[test]
fn testing_mgi_full_header() {
    for ds_itr_tmp in 1..5{
        
        let ds_itr_tmp = if ds_itr_tmp == 3{
            8
        }else if ds_itr_tmp == 4 {11}
        else{ds_itr_tmp};

        let ds_itr_in = if ds_itr_tmp == 11{
            1
        }else{ds_itr_tmp};
        
        let ds_itr_fc = if ds_itr_tmp == 11{
            1
        }else{ds_itr_tmp};
        
        let ds_itr_ex = ds_itr_tmp;
        
        let read1_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_1.fq.gz", ds_itr_in, ds_itr_fc));
        let read2_file_path : String = String::from(format!("testing_data/input/ds0{}/L01/FC0{}_L01_read_2.fq.gz", ds_itr_in, ds_itr_fc));
        let ext = "tsv";
        let sample_sheet_file_path : String = String::from(format!("testing_data/expected/ds0{}/sample_sheet_expected.{ext}", ds_itr_ex));
        let lane = String::from("L01");
        let comprehensive_scan = true;
        
        for allowed_mismatches in 0..1 {
            
            let ouput_dir = format!("testing_data/output/full_header/out_real-{}/", ds_itr_in);
            if PathBuf::from(&ouput_dir).exists() {
                fs::remove_dir_all(&ouput_dir).unwrap();
            }
            let original_path = format!("testing_data/expected/full_mgi_header-{}/", ds_itr_tmp);
        
            let command = "target/debug/mgikit";
            let mut my_args: Vec<String> = vec!["demultiplex".to_string(),
                                                "-f".to_string(),
                                                read1_file_path.to_string(), 
                                                "-r".to_string(), 
                                                read2_file_path.to_string(), 
                                                "-s".to_string(), 
                                                sample_sheet_file_path.to_string(), 
                                                "--lane".to_string(), 
                                                lane.to_string(), 
                                                "--writing-buffer-size".to_string(), 
                                                "131072".to_string(), 
                                                "-o".to_string(),
                                                ouput_dir.to_string(), 
                                                "-m".to_string(), 
                                                format!("{}", allowed_mismatches), 
                                                "--force".to_string(),
                                                "-t".to_string(),
                                                "1".to_string()];
                        
            if comprehensive_scan{
                my_args.push("--comprehensive-scan".to_string());
            }
            
            my_args.push("--disable-illumina".to_string());
            my_args.push("--mgi-full-header".to_string());
            
            if ds_itr_tmp <  11{
                my_args.push("--all-index-error".to_string());                
            }
            my_args.push("--validate".to_string());
            println!("{:?}", &my_args);


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
      
}

struct TestCleanup;

impl Drop for TestCleanup {
    fn drop(&mut self) {
        for entry in fs::read_dir("testing_data/output/").unwrap() {
            let entry = entry.unwrap();
            let file_type = entry.file_type().unwrap();
    
            let entry_path = entry.path();
            if file_type.is_dir() {
                // Delete the directory itself
                fs::remove_dir_all(entry_path).unwrap();
            } else {
                // Delete the file
                fs::remove_file(entry_path).unwrap();
            }
        }
    }
}