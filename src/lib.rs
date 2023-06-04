
use std::collections::{HashMap};
use std::path::Path;
use std::fs::File;
use std::time::Instant;
use flate2::read::{MultiGzDecoder};
use std::io::{self, BufRead};
use std::fs::OpenOptions;
use std::io::BufReader;
use std::io::Write;
use chrono::prelude::*;
use std::fs;
use gzp::{
    deflate::Mgzip,
    par::compress::{ParCompress, ParCompressBuilder},
    ZWriter,
};



// my modules
mod variables;
use crate::variables::*;
mod samplesheet_utils;
use crate::samplesheet_utils::*;
mod sequence_utils;
use crate::sequence_utils::*;
mod index_dic;
use crate::index_dic::*;



fn write_reads(reads_buffer: &Vec<String>, output_file: &String) {
    //return;
    //println!("Writing {}", output_file);
    //let append = true;
    let output_file_path = Path::new(output_file);
    //let chunksize = 64 * (1 << 10) * 2;
    let outfile = OpenOptions::new()
        .append(true)
        .create(true)
        .open(output_file_path)
        .expect("couldn't create output");

    let mut writer: ParCompress<Mgzip> = ParCompressBuilder::new().from_writer(outfile);

    for read_inf in reads_buffer {
        writer.write_all(&read_inf.as_bytes()).unwrap();
    }
    writer.finish().unwrap();
}

fn check_file(path_str: &String) {
    if !Path::new(path_str).exists() {
        panic!("File is not accessible: {}", path_str);
    }
}

fn create_folder(path_str: &String) {
    println!("A directory is created: {}", path_str);
    fs::create_dir(path_str).unwrap();
}


fn convert_read_header_to_illumina(mgi_read_header: &String, umi: &String) -> String {
    
    //println!("header: {}", &mgi_read_header);
    let mut header_itr = 0;
    let mut l_position = 0;
    for header_chr in mgi_read_header.chars().rev() {
        if header_chr == 'L' {
            l_position = header_itr;
            break;
        }

        header_itr += 1;
    }

    let mut output = String::from(":");
    output.push_str(&mgi_read_header[1..(mgi_read_header.len() - l_position - 1)]);
    //output.push_str(&"V350015219");
    
    output.push(':');
    output.push_str(
        &mgi_read_header
            [(mgi_read_header.len() - l_position)..(mgi_read_header.len() - l_position + 1)],
    );
    output.push(':');
    output.push_str(
        &mgi_read_header[(mgi_read_header.len() - l_position + 9)..(mgi_read_header.len() - 3)]
            .parse::<i32>()
            .unwrap()
            .to_string(),
    );
    output.push(':');
    output.push_str(
        &mgi_read_header
            [(mgi_read_header.len() - l_position + 2)..(mgi_read_header.len() - l_position + 5)]
            .parse::<i32>()
            .unwrap()
            .to_string(),
    );
    output.push(':');
    output.push_str(
        &mgi_read_header
            [(mgi_read_header.len() - l_position + 6)..(mgi_read_header.len() - l_position + 9)]
            .parse::<i32>()
            .unwrap()
            .to_string(),
    );
    output.push_str(&umi);
    output.push(' ');
    output.push_str(&mgi_read_header[(mgi_read_header.len() - 2)..(mgi_read_header.len() - 1)]);
    output.push_str(&":N:0:");
    output
}


pub fn write_general_info_report(sample_information:&Vec<Vec<String>>, sample_statistics:&Vec<Vec<u64>>, kept_samples: &Vec<usize>, run:&String, lane:&String, output_file:String){
     //Mb Total Yield: total bases as cnt of r1 + r2 	
    //M Total Clusters: number of reads	
    //% bases ≥ Q30	r1 and r2 bases that are greater than 30
    //Mean Quality:	r1qc + r2qc / cnt or bases in r1 and r2
    //% Perfect Index: index with 0 mismatches / all reads.

    /*
            sample_statistics[sample_id]:
            0: r1 count of qs > 30
            1: r2 count os qs > 30
            2: r3 count of (barcode) with qs > 30

            3: r1 count of bases
            4: r2 count of bases
            5: r3 count of bases

            6: qs of r1
            7: qs of r2
            8:  qs of r3
            9: read count
            from 10 to whatever mismatches allowed: is the number oif reads with the iteration mismatch
            
        */
    let mut out_str = String::from("#sample general info\nSample ID\tM Clusters\tMb Yield ≥ Q30\t% R1 Yield ≥ Q30\t% R2 Yield ≥ Q30\t% R3 Yield ≥ Q30\t% Perfect Index\n");
    let non_sample = vec!["Undetermined".to_string(), "Ambiguous".to_string()];
    let mut all_q_30;
    let exclude_non_sample = true;
    let mut outfile;
    let mut curr_sample_stat: Vec<f64>;
    let mut tmp_val;
    let mut sample_id;
    let mut lane_statistics:  Vec<f64> = vec![0 as f64, 0 as f64, 0 as f64, 0 as f64, 0 as f64];
    let mut sample_statistics_copy = sample_statistics.clone();
    for sample_id_itr in 0..sample_statistics_copy.len() {
        sample_statistics_copy[sample_id_itr].push(sample_id_itr as u64);
    }
    
    sample_statistics_copy.sort_by(|a, b| b[9].cmp(&a[9]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap())));
    for sample_id_itr in 0..sample_statistics_copy.len(){
        sample_id = sample_statistics_copy[sample_id_itr].last().unwrap().clone() as usize;
        if kept_samples.len() > 0 && !kept_samples.contains(&sample_id){
            continue;
        }

        curr_sample_stat = sample_statistics_copy[sample_id_itr].iter().map(|x| x.clone() as f64).collect();

        out_str.push_str(&sample_information[sample_id][SAMPLE_COLUMN]);
        out_str.push('\t');

        // M clusters
        out_str.push_str(&format!("{}", curr_sample_stat[9] / 1000000.0));
        out_str.push('\t');
        all_q_30 = curr_sample_stat[0] + curr_sample_stat[1];
        lane_statistics[0] += curr_sample_stat[3] + curr_sample_stat[4];
        lane_statistics[1] += curr_sample_stat[9];
        lane_statistics[2] += all_q_30;
        
        // Mb Yiled > 30
        out_str.push_str(&format!("{}", all_q_30 / 1000000.0));
        out_str.push('\t');
        
        // % R1 > 30
        if curr_sample_stat[3] == 0.0{
            tmp_val = 0.0;
        }
        else{
            tmp_val = curr_sample_stat[0] / curr_sample_stat[3];
        }
        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));
        out_str.push('\t');
            
        // % R2 > 30
        if curr_sample_stat[4] == 0.0{
            tmp_val = 0.0;
        }
        else{
            tmp_val = curr_sample_stat[1] / curr_sample_stat[4];
        }

        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));
        
        // % R3 > 30
        out_str.push('\t');

        if curr_sample_stat[5] == 0.0{
            tmp_val = 0.0;
        }
        else{
            tmp_val = curr_sample_stat[2] / curr_sample_stat[5];
        }
        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));
        out_str.push('\t');
        
        
        // Perfect index
        if curr_sample_stat[9] == 0.0{
            tmp_val = 0.0;
        }
        else{
            tmp_val = curr_sample_stat[10] / curr_sample_stat[9];
        }
        
        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));

        lane_statistics[3] += curr_sample_stat[6] + curr_sample_stat[7];
        
        if !exclude_non_sample || !non_sample.contains(&sample_information[sample_id][SAMPLE_COLUMN].split("_").map(|x| x.to_string()).collect::<Vec<String>>()[0].to_lowercase()){
            lane_statistics[4] += curr_sample_stat[10];
        }
        
        out_str.push('\n');       

    }
    
    //output_file_path = Path::new(&output_file);
    let mut final_out_str = String::from("#Lane statistics\nRun ID-Lane\tMb Total Yield\tM Total Clusters\t% bases ≥ Q30\tMean Quality\t% Perfect Index\n");
    final_out_str.push_str(&run);
    final_out_str.push('-');
    final_out_str.push_str(&lane);
    final_out_str.push('\t');

    final_out_str.push_str(&format!("{}", lane_statistics[0] / 1000000.0));
    final_out_str.push('\t');

    final_out_str.push_str(&format!("{}", lane_statistics[1] / 1000000.0));
    final_out_str.push('\t');


    final_out_str.push_str(&format!("{:.3?}", lane_statistics[2] / lane_statistics[0] * 100.0));
    final_out_str.push('\t');
    final_out_str.push_str(&format!("{:.3?}", lane_statistics[3] / lane_statistics[0] ));
    final_out_str.push('\t');
    final_out_str.push_str(&format!("{:.3?}", lane_statistics[4] / lane_statistics[1] * 100.0));
    final_out_str.push('\n');
    
    final_out_str.push_str(&out_str);

    outfile = File::create(Path::new(&output_file)).expect("couldn't create output");
    outfile.write_all(&final_out_str.as_bytes()).unwrap();
   

}



pub fn write_index_info_report(sample_information:&Vec<Vec<String>>, sample_mismatches:&Vec<Vec<u64>>, kept_samples: &Vec<usize>, max_mismatches:usize, output_file:String){
    let mut report_str = String::from("sample");

    for i in 0..max_mismatches {
        report_str.push('\t');
        report_str.push_str(&i.to_string());
        report_str.push_str(&"-mismatches");
    }
    report_str.push('\n');
    let mut sample_mismatches_copy = sample_mismatches.clone();
    for sample_id_itr in 0..sample_mismatches_copy.len() {
        sample_mismatches_copy[sample_id_itr].push(sample_id_itr as u64);
    }
    
    //sample_mismatches_copy.sort_by(|a, b| b[0].cmp(&a[0]));
    sample_mismatches_copy.sort_by(|a, b| b[0].cmp(&a[0]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap())));
    //sample_mismatches_copy.sort_by(|a, b| b[0].cmp(&a[0]));
    for sample_id_itr in 0..sample_mismatches_copy.len() {
        let sample_index = *sample_mismatches_copy[sample_id_itr].last().unwrap() as usize;
        if kept_samples.len() > 0 && ! kept_samples.contains(&sample_index) {
            continue;
        } 
        report_str.push_str(&sample_information[sample_index][SAMPLE_COLUMN]);
        for cnt in 1..(max_mismatches + 1) {
            if cnt < sample_mismatches_copy[sample_id_itr].len() {
                report_str.push('\t');
                report_str.push_str(&sample_mismatches_copy[sample_id_itr][cnt].to_string());
            } else {
                report_str.push_str("\t0");
            }
        }
        report_str.push('\n');
    }

    //println!("report: {}", output_file_path.display());
    let mut outfile = File::create(Path::new(&output_file)).expect("couldn't create output");
    outfile.write_all(&report_str.as_bytes()).unwrap();
        
}


/// demultiplex docs
/// This is the main funciton in this crate, which demultiplex fastq single/paired end files and output samples' fastq and quality and run statistics reports.
/// 
/// # Arguments
///
/// * `input_folder_path` - A string that holds the path to teh directory oif the run generated by MGI machines (this usually referes to the Lane directory).
/// * `read1_file_path` - A string that holds the path to the farward fastq file generated by the sequncing machine to be demultiplexed.
/// * `read_file_path` - A string that holds the path to the reveresed fastq file generated by the sequncing machine to be demultiplexed.
///
/// # Examples

pub fn demultiplex(
    input_folder_path: &String,
    read2_file_path: &String,
    read1_file_path: &String,
    sample_sheet_file_path: &String,
    ouput_dir: &String,
    report_dir: &String,
    allowed_mismatches: usize,
    template: &String,
    arg_i7_rc: bool,
    arg_i5_rc: bool,
    arg_lane: &String,
    arg_instrument: &String,
    arg_run: &String,
    disable_illumina_format: bool,
    //single_read_input: bool,
    keep_barcode: bool,
    writing_threshold: usize,
    read_merging_threshold: usize,
    comprehensive_scan: bool,
    undetermined_label: &String,
    ambiguous_label: &String,
    force: bool,
    report_limit: usize

) {
    // Validate input data
    let mut single_read_input = false;
    let read1_file_name_suf = "_read_1.fq.gz";
    let read2_file_name_suf = "_read_2.fq.gz";
    let info_file = "BioInfo.csv";
    
    let mut info_file_path = String::new();
    let mut paired_read_file_path_final  = read1_file_path.clone();
    let mut read_barcode_file_path_final = read2_file_path.clone();
    
    let mut instrument = String::new();
    let mut run = String::new();
    let mut lane = String::new();
    if input_folder_path.len() > 0{
        let entries = fs::read_dir(input_folder_path).unwrap()
        .map(|res| res.map(|e| e.path()))
        .collect::<Result<Vec<_>, io::Error>>().unwrap();
        for path in entries{
            if path.file_name().unwrap().to_str().unwrap().ends_with(read1_file_name_suf){
                let tmp : Vec<&str> = path.file_name().unwrap().to_str().unwrap().split("_").collect();
                if tmp.len() > 1{
                    lane = tmp[1].to_string().trim().to_string();
                    if ! lane.starts_with("L0"){
                        lane = String::new();
                    }
                }
                paired_read_file_path_final = path.to_str().unwrap().to_string();
            }else if path.file_name().unwrap().to_str().unwrap().ends_with(read2_file_name_suf){
                read_barcode_file_path_final = path.to_str().unwrap().to_string();
            }else if path.file_name().unwrap().to_str().unwrap().ends_with(info_file){
                info_file_path = path.to_str().unwrap().to_string();
            }
        }
        
        println!("R1: {}", paired_read_file_path_final);
        println!("R2: {}", read_barcode_file_path_final);
        println!("info: {}", info_file_path);
    }

    if info_file_path.len() > 0{
        let file = File::open(info_file_path).unwrap();
        for line in io::BufReader::new(file).lines() {
                if let Ok(inf) = line {
                    //println!("{}", inf);
                    if inf.starts_with("Machine ID,"){
                        let tmp : Vec<&str> = inf.split(",").collect();
                        instrument = tmp[1].to_string().trim().to_string();
                    }else if inf.starts_with("Sequence Date") || inf.starts_with("Sequence Start Date"){
                        let tmp : Vec<&str> = inf.split(",").collect();
                        run = tmp[1].to_string().trim().replace("-", "");
                    }else if inf.starts_with("Sequence Time") || inf.starts_with("Sequence Start Time"){
                        let tmp : Vec<&str> = inf.split(",").collect();
                        run.push_str(&tmp[1].to_string().trim().replace(":", ""));
                    }
                }
            
        }
    }

    if arg_instrument.len() > 0{
        instrument = arg_instrument.clone();
    }
    if arg_run.len() > 0{
        run = arg_run.clone();
    }
    if arg_lane.len() > 0{
        lane = arg_lane.clone();
    }
       

    if read_barcode_file_path_final.len() == 0{
        //Single End reads
        read_barcode_file_path_final = paired_read_file_path_final;
        paired_read_file_path_final = String::new();
        single_read_input = true;
        println!("Single end read input was detected!");
        //if input_barcode_length == 0{
        //    panic!("Barcode length should be greater than 0 when processing single end reads");
        //}
    }else{
        println!("Paired ended read input was detected!");
    }
    

    


    if read_barcode_file_path_final.len() == 0 {
        panic!("Input reads are invalid! check the path {}", read_barcode_file_path_final);
    }

    check_file(&read_barcode_file_path_final);
    
    if !single_read_input {
        if paired_read_file_path_final.len() == 0 {
            panic!("Input reads are invalid! check the path {}", paired_read_file_path_final);
        }
        check_file(&paired_read_file_path_final);
    }

    if sample_sheet_file_path.len() == 0 {
        panic!("Sample sheet file is invalid!");
    }
    check_file(sample_sheet_file_path);



    let mut output_directory = ouput_dir.clone();
    if output_directory.len() == 0 {
        let local: DateTime<Local> = Local::now();
        output_directory = String::from("mgiKit_");
        output_directory.push_str(&local.format("%Y%m%dT%H%M%S").to_string());
        output_directory.push('/');
    }else{
        if output_directory.chars().last().unwrap() != '/'{
            output_directory.push('/');
        }
    }

    let mut report_dir_local = report_dir.clone();
    if report_dir_local.len() == 0 {
        println!("The same output directory will be used for reports.");
        report_dir_local = output_directory.clone();
    }
    if report_dir_local.chars().last().unwrap() != '/'{
        report_dir_local.push('/');
    }

    if (Path::new(&output_directory).is_dir() || Path::new(&report_dir_local).is_dir()) && !force {
        panic!(
                "One or both of output directly and report directly exists. Use --froce to overwrite their data: {} and {}",
                output_directory, report_dir_local
            );
        
    } 
    
    if Path::new(&output_directory).is_dir(){
        println!(
            "Output directory exists. Data will be overwritten at: {}.",
            output_directory
        );
    }else {
        create_folder(&output_directory);
    }

    if Path::new(&report_dir_local).is_dir(){
        if report_dir_local != output_directory{
            println!(
                "Report directory exists. Data will be overwritten at: {}.",
                report_dir_local
            );
        }
        
    }else {
        create_folder(&report_dir_local);
    }

    /*if input_barcode_length == 0 {
        println!("barcode length will be calcualated as the length's difference between R2 and R1!");
    }else{
        println!("barcode length is provided by the user as ({} bp)", input_barcode_length);
    }*/
    
    
    println!("Output directory: {}", output_directory);
    println!("Reports directory: {}", report_dir_local);
    println!("Instrumnet: {}", instrument);
    println!("Run: {}", run);
    println!("Lane: {}", lane);
    println!("Comprehensive scan mood: {}", comprehensive_scan);
   

    let mut i7_rc = arg_i7_rc;
    let i5_rc = arg_i5_rc;

    let trim_barcode = !keep_barcode;
    let illumina_format = !disable_illumina_format;

    if template.len() > 0 {
        println!(
            "General template is provided and will be used for all samples: {}",
            template
        );
        if i7_rc {
            println!("i7 will be converted to the reverse complement!");
        }

        if i5_rc {
            println!("i5 will be converted to the reverse complement!");
        }
    } else {
        println!("Template will be used from sample/index map file.");
        i7_rc = false;
        //i5_rc = false;
    }

    println!(
        "Allowed mismatches when finding the index are: {}",
        allowed_mismatches
    );

    println!(
        "Reads that match with multiple samples will be saved in ambiguous_read1/2.fastq file"
    );

    println!("The paired read files are assumed to contain the paired reads in the same order!");

    let mut illumina_header_template_p1 = String::new();
    let mut output_file_format_r2;
    let mut output_file_format_r1 = String::new();

    if illumina_format {
        if lane.len() == 0 || instrument.len() == 0 || run.len() == 0 {
            panic!("Instrument id and run number are required for QC reports and when Illumina format is requested!")
        }
        
        illumina_header_template_p1 = "@".to_string();
        illumina_header_template_p1.push_str(&instrument);
        illumina_header_template_p1.push(':');
        illumina_header_template_p1.push_str(&run);

        output_file_format_r2 = lane.clone();
        
            
        if single_read_input{
            output_file_format_r2.push_str(&"_R1_001.fastq.gz");

        }else{
            output_file_format_r1 = lane.clone();
            output_file_format_r1.push_str(&"_R1_001.fastq.gz");
            output_file_format_r2.push_str(&"_R2_001.fastq.gz");
        }

        println!("Read header and Output files: Illumina format.");

    } else {
        println!("Read header and Output files: MGI format.");
        if lane.len() == 0 {
            panic!("Lane number is required for QC reports!")
        }

        if single_read_input{
            output_file_format_r2 = "_R1.fastq.gz".to_string();
        }else{
            output_file_format_r1 = "_R1.fastq.gz".to_string();
            output_file_format_r2 = "_R2.fastq.gz".to_string();
        }
        //extract_umi = false;
    }

    println!("Trim Barcode: {}", trim_barcode);

    // writing_threshold: usize, read_merging_threshold
    println!(
        "Writing configuration: {} -  {}",
        writing_threshold, read_merging_threshold
    );

    //let mut sample_output_buffer_read_barcode: HashMap<String, (Vec<u32>, String)> = HashMap::new();
    //let mut sample_output_buffer_paired_read: HashMap<String, (Vec<u32>, String)> = HashMap::new();
    // parse sample/index file and get all mismatches
    
    let tmp_res = parse_sample_index(Path::new(sample_sheet_file_path), &template, i7_rc, i5_rc).unwrap();
    let all_template_data = tmp_res.0;
    let mut sample_information = tmp_res.1;
    let project_samples = tmp_res.2;
    let undetermined_label_id = sample_information.len();
    let ambiguous_label_id = sample_information.len() + 1;
    sample_information.push(vec![undetermined_label.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), String::from(".")]);
    sample_information.push(vec![ambiguous_label.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), String::from(".")]);

    let mut mismatches_dic_i7: Vec<HashMap<String, (Vec<&String>, usize)>>  = Vec::new();
    let mut mismatches_dic_i5: Vec<HashMap<String, (Vec<&String>, usize)>>  = Vec::new();
    
    for template_details in &all_template_data {
        mismatches_dic_i7.push(get_all_mismatches(&template_details.1, allowed_mismatches));
        mismatches_dic_i5.push(get_all_mismatches(&template_details.2, allowed_mismatches));
    }
        

    //println!("keys: {}", template_index_sample_information[0].keys().len());

    //println!("{:?}", all_i7s);
    //println!("{:?}", all_i5s);
    //return;

    //let mut read_barcode_data = BGZFReader::new(File::open(path_to_read2).expect("Could not open the file"));
    //let mut read_barcode_data = BufReader::new(ConsecutiveReader:: from_path(path_to_read2, 0).expect("Could not open the file"));
    // both above works fine on Bgzip
    //let mut read_barcode_data = BufReader::new(GzDecoder::new(File::open(path_to_read2).expect("Could not open the file")));
    // this above works for gzip files.

    // this above works for gzip files.
    /*
    let mut lines = 0;
    loop
    {
        let mut line = String::new();
        let read_bytes = read_barcode_data.read_line(&mut line).unwrap();
        println!("{} - {} ", lines, read_bytes);
        if read_bytes == 0 {
            println!("bytes: {}", read_bytes);
            break;
        }
        lines += 1;
    }

    return;
    */

    let mut tmp_str = String::new();
    let barcode_read_length: usize;
    let paired_read_length;
    let mut read2_has_sequence = true;
    let mut read_barcode_data = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path_final)).expect("Could not open the file")
    ));
    read_barcode_data.read_line(&mut tmp_str).unwrap();
    tmp_str = String::new();
    read_barcode_data.read_line(&mut tmp_str).unwrap();
    barcode_read_length = tmp_str.chars().count() - 1;

    let mut paired_read_data = if !single_read_input {
        Some(BufReader::new(MultiGzDecoder::new(
            File::open(Path::new(&paired_read_file_path_final)).expect("Could not open the file")
        )))
    }else{
        None
    };
    
    match paired_read_data.as_mut() {
        Some(paired_read_data_buff) => {
            paired_read_data_buff.read_line(&mut tmp_str).unwrap();
            tmp_str = String::new();
            paired_read_data_buff.read_line(&mut tmp_str).unwrap();
            paired_read_length = tmp_str.chars().count() - 1;
        },
        None => {
            paired_read_length = 0;
        }
    };
    
    println!("The length of the read with barcode is: {}", barcode_read_length);
    println!("The length of the paired read is: {}", paired_read_length);
    if paired_read_length >= barcode_read_length{
        read2_has_sequence = false;
        println!("It is assumed that read2 contains barcode only without read sequence!");    
    }

    read_barcode_data = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path_final)).expect("Could not open the file")
    ));

    paired_read_data = if !single_read_input {
        Some(BufReader::new(MultiGzDecoder::new(
            File::open(Path::new(&paired_read_file_path_final)).expect("Could not open the file")
        )))
    }else{
        None
    };

    let mut read_barcode_header;
    let mut read_barcode_seq;
    let mut read_barcode_plus;
    let mut read_barcode_quality_score;

    let mut paired_read_header;
    let mut paired_read_seq;
    let mut paired_read_plus;
    let mut paired_read_quality_score;

    let mut extract_umi;
    let mut read_bytes;

    let mut out_read_barcode_buffer: Vec<Vec<String>> = Vec::new();
    let mut out_paired_read_buffer: Vec<Vec<String>> = Vec::new();
    let mut sample_mismatches: Vec<Vec<u64>> = Vec::new();
    let mut sample_output_path: Vec<(String, String)> = Vec::new();
    let mut sample_statistics: Vec<Vec<u64>> = Vec::new();
    
    for i in 0..sample_information.len(){
        out_read_barcode_buffer.push(vec![String::new()]);
        out_paired_read_buffer.push(vec![String::new()]);
        sample_mismatches.push(vec![0; 2 * allowed_mismatches + 2]);
        
        
        let mut tmp = String::from(&output_directory);
        tmp.push_str(&sample_information[i][SAMPLE_COLUMN]);
        
        if i < undetermined_label_id{
            if illumina_format{
                tmp.push_str(&"_S".to_string());
                tmp.push_str(&(i + 1).to_string());
                tmp.push('_');
    
            }      
            tmp.push_str(&output_file_format_r2);
        }else{
            if illumina_format{
                tmp.push('_');
            }
            tmp.push_str(&output_file_format_r2);
        }

        if Path::new(&tmp).exists(){
            fs::remove_file(&tmp).unwrap();
        }
        
        let mut tmp1 = String::from(&output_directory);
        tmp1.push_str(&sample_information[i][SAMPLE_COLUMN]);
        
        if i < undetermined_label_id{
            if illumina_format{
                tmp1.push_str(&"_S".to_string());
                tmp1.push_str(&(i + 1).to_string());
                tmp1.push('_');
            }
            tmp1.push_str(&output_file_format_r1);    
        }else{
            if illumina_format{
                tmp1.push('_');
            }
            tmp1.push_str(&output_file_format_r1);
        }
        
        if Path::new(&tmp1).exists(){
            fs::remove_file(&tmp1).unwrap();
        }
        sample_output_path.push((tmp, tmp1));
        sample_statistics.push(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

    }

    //let mut curr_template;
    let mut sample_id;
    let mut sample_info;
    let mut barcode_length=0;
    
    let mut indexes_info;
    let mut i7_read;
    let mut i5_read = String::new();
    let mut read_cntr:u64 = 0;
    let mut curr_mismatch;
    let mut curr_umi = String::new();
    let mut curr_barcode = String::new();
    let mut matching_samples: Vec<(String, u32)> = Vec::new();
    let mut undetermined_barcodes: HashMap<String, u32> = HashMap::new();
    let mut ambiguous_barcodes: HashMap<String, u32> = HashMap::new();
    //println!("start looping");
    let mut last_piece;
    let mut template_itr;

    let mut flowcell = String::new(); 
    //let mut flowcell = String::from("V350015219");
    
    let mut all_i7_mismatches;
    let mut all_i5_mismatches;
    let mut tmp_cntr;
    let shift = if single_read_input{0}else{1};

    let start = Instant::now();
    let dur;
    //print_type_of(&read_cntr);
    loop {
        //println!("**************");
        /*if  read_cntr % 5000000 == 0 {
            dur = start.elapsed();
            println!("{} reads  {} seconds", read_cntr, dur.as_secs());
            //break;
        }*/
        //println!("{} reads", lines);
        paired_read_header = String::new();
        paired_read_seq = String::new();
        paired_read_plus = String::new();
        paired_read_quality_score = String::new();
        read_barcode_header = String::new();
        read_barcode_seq = String::new();
        read_barcode_plus = String::new();
        read_barcode_quality_score = String::new();
        curr_mismatch = 0;

        //let mut read_barcode_info = String::new();
        //let start = Instant::now();
        
        read_bytes = read_barcode_data.read_line(&mut read_barcode_header).unwrap();
        //println!("header: {} ", read_barcode_header);
        if read_bytes == 0 {
            //println!("bytes: {}", read_bytes);
            break;
        }

        if flowcell.len() == 0{
            let v: Vec<&str> = read_barcode_header.split('L').collect();
            flowcell = v[0][1..].to_string();
        }

        read_barcode_data.read_line(&mut read_barcode_seq).unwrap();
        //println!("read: *{}*", read_barcode_seq);
        //println!("read length: {}  ->  {}", read_barcode_seq.len(), read_barcode_seq.chars().count());

        read_barcode_data.read_line(&mut read_barcode_plus).unwrap();
        read_barcode_data.read_line(&mut read_barcode_quality_score).unwrap();
        //let dur = start.elapsed();
        //t_r2_read += dur.as_nanos();
        //println!("{}", t_r2_read);
        //continue;

        sample_id = sample_information.len();
        matching_samples.clear();
        template_itr = 0;

        //let start = Instant::now();
        
        
        {
            //println!("-------------");
            for template_details in &all_template_data {
                sample_info = &template_details.4;
                indexes_info = template_details.6;
                //println!("index: {:?}", indexes_info);
                all_i7_mismatches = &mismatches_dic_i7[template_itr];
                all_i5_mismatches = &mismatches_dic_i5[template_itr];
                


                if indexes_info[6] == 1 {
                    extract_umi = true;
                } else {
                    extract_umi = false;
                }

                i7_read = read_barcode_seq[(read_barcode_seq.len() - indexes_info[2] - 1)
                    ..(read_barcode_seq.len() - indexes_info[2] + indexes_info[1] - 1)]
                    .to_string();
                
                //println!("-----------------");
                
                //println!("template: {:?}", indexes_info[1]);
                //println!("bytes: {}", read_bytes);

                curr_barcode = i7_read.clone();
                if template_details.5 {
                    i5_read = read_barcode_seq[(read_barcode_seq.len() - indexes_info[5] - 1)
                            ..(read_barcode_seq.len() - indexes_info[5] + indexes_info[4] - 1)]
                            .to_string();
                    curr_barcode.push('+');
                    curr_barcode.push_str(&i5_read);
                }
                
                
                match all_i7_mismatches.get(&i7_read) {
                    Some(i7_matches) => {
                        if template_details.5 {
                            match all_i5_mismatches.get(&i5_read) {
                                Some(i5_matches) => {
                                    //println!("{:?} and {:?}", i7_matches, i5_matches);
                                    for i7_match_itr in 0..i7_matches.0.len() {
                                        for i5_match_itr in 0..i5_matches.0.len() {
                                            curr_mismatch = i7_matches.1 + i5_matches.1;
                                            if curr_mismatch <= allowed_mismatches{
                                                //println!("{}", i7_matches.0[i7_match_itr]);
                                                match sample_info
                                                .get(i7_matches.0[i7_match_itr])
                                                .unwrap()
                                                .1
                                                .get(i5_matches.0[i5_match_itr])
                                                {
                                                    Some(i5_info) => {
                                                        if sample_id >= sample_information.len() {
                                                            sample_id = *i5_info;
                                                            barcode_length = indexes_info[9];
                                                            if illumina_format {
                                                                if extract_umi {
                                                                    curr_umi = String::from(":");
                                                                    if i7_rc {
                                                                        curr_umi.push_str(
                                                                            &reverse_complement(
                                                                                &read_barcode_seq[(read_barcode_seq
                                                                                    .len()
                                                                                    - indexes_info[8]
                                                                                    - 1)
                                                                                    ..(read_barcode_seq.len()
                                                                                        - indexes_info
                                                                                            [8]
                                                                                        + indexes_info
                                                                                            [7]
                                                                                        - 1)]
                                                                                    .to_string(),
                                                                            )
                                                                            .unwrap(),
                                                                        );
                                                                    } else {
                                                                        curr_umi.push_str(
                                                                            &read_barcode_seq[(read_barcode_seq.len()
                                                                                - indexes_info[8]
                                                                                - 1)
                                                                                ..(read_barcode_seq.len()
                                                                                    - indexes_info[8]
                                                                                    + indexes_info[7]
                                                                                    - 1)],
                                                                        );
                                                                    }
                                                                }
                                                            }
                                                        } else {
                                                            sample_id = ambiguous_label_id;
                                                            break;
                                                        }
                                                    }
                                                    None => {
                                                        //break
                                                    },
                                                };
                                            }
                                        }
                                    }
                                    /*if sample_id >= sample_information.len() {
                                        sample_id = undetermined_label_id;
                                    }*/
                                },
                                None => {
                                    /*if sample_id >= sample_information.len() {
                                        sample_id = undetermined_label_id;
                                    }*/
                                }
                            }
                        } else {
                            if i7_matches.0.len() > 1 || sample_id < sample_information.len() - 2 {
                                sample_id = ambiguous_label_id;
                            } else {
                                sample_id = template_details.4.get(i7_matches.0[0]).unwrap().0;
                                barcode_length = indexes_info[9];
                                curr_mismatch = i7_matches.1;
                                if illumina_format {
                                    curr_barcode = i7_read.clone();
                                    if extract_umi {
                                        curr_umi = String::from(":");
                                        curr_umi.push_str(
                                            &read_barcode_seq[(read_barcode_seq.len() - indexes_info[8] - 1)
                                                ..(read_barcode_seq.len() - indexes_info[8]
                                                    + indexes_info[7]
                                                    - 1)],
                                        );
                                        if i7_rc {
                                            curr_umi.push_str(
                                                &reverse_complement(
                                                    &read_barcode_seq[(read_barcode_seq.len()
                                                        - indexes_info[8]
                                                        - 1)
                                                        ..(read_barcode_seq.len() - indexes_info[8]
                                                            + indexes_info[7]
                                                            - 1)]
                                                        .to_string(),
                                                )
                                                .unwrap(),
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    },
                    None => {
                        /*if sample_id >= sample_information.len() {
                            sample_id = undetermined_label_id;
                        }*/
                    }
                };
                
                //println!("{} - {} ->  {:?}",undetermined_label_id, sample_id, indexes_info);
                
                if (sample_id != undetermined_label_id && sample_id < sample_information.len()) && !comprehensive_scan {
                    break;
                }

                
                template_itr += 1;
            }
        }
       
        if sample_id >= sample_information.len() {
            sample_id = undetermined_label_id;
            barcode_length = all_template_data[0].6[9];   
        }

        if sample_id == ambiguous_label_id {
            curr_mismatch = 0;
            match ambiguous_barcodes.get_mut(&curr_barcode) {
                Some(barcode_info) => {
                    *barcode_info += 1;
                }
                None => {
                    ambiguous_barcodes.insert(curr_barcode.clone(), 1);
                }
            };
            
        }

        //if illumina_format && sample_id != undetermined_label_id && sample_id != ambiguous_label_id {
        
        if read2_has_sequence && illumina_format {
                let illumina_headers_elements =
                convert_read_header_to_illumina(&read_barcode_header, &curr_umi);
            read_barcode_header = illumina_header_template_p1.clone();
            read_barcode_header.push_str(&illumina_headers_elements);
            read_barcode_header.push_str(&curr_barcode);
            read_barcode_header.push('\n');
        }

        if sample_id == undetermined_label_id {
            curr_mismatch = 0;
            match undetermined_barcodes.get_mut(&curr_barcode) {
                Some(barcode_info) => {
                    *barcode_info += 1;
                }
                None => {
                    undetermined_barcodes.insert(curr_barcode.clone(), 1);
                }
            };
        }


        
        match paired_read_data.as_mut() {
            Some(paired_read_data_buff) => {
                paired_read_data_buff.read_line(&mut paired_read_header).unwrap();
                paired_read_data_buff.read_line(&mut paired_read_seq).unwrap();
                paired_read_data_buff.read_line(&mut paired_read_plus).unwrap();
                paired_read_data_buff.read_line(&mut paired_read_quality_score).unwrap();

                //if barcode_length == 0 {
                //    barcode_length = read_barcode_seq.len() - paired_read_seq.len();
                //    println!("Detected barcode length is {} bp (R2: {}, R1: {})", barcode_length, read_barcode_seq.len() - 1, paired_read_seq.len() - 1)
                    //println!("ZZZ: {} {} {}", read_barcode_seq.len(), paired_read_seq.len(), barcode_length);
                    //println!("V: {:}", flowcell);
                //}
                

                //if illumina_format && sample_id != undetermined_label_id && sample_id != ambiguous_label_id {
                if illumina_format {
                    let illumina_headers_elements =
                        convert_read_header_to_illumina(&paired_read_header, &curr_umi);
                    paired_read_header = illumina_header_template_p1.clone();
                    paired_read_header.push_str(&illumina_headers_elements);
                    paired_read_header.push_str(&curr_barcode);
                    paired_read_header.push('\n');
                }
            },
            None    => {},
        }
        
        
        

        /*
            sample_statistics[sample_id]:
            0: r1 count of qs > 30
            1: r2 count os qs > 30
            2: r3 count of (barcode) with qs > 30

            3: r1 count of bases
            4: r2 count of bases
            5: r3 count of bases

            6: qs of r1
            7: qs of r2
            8:  qs of r3
            9: read count
            from 10 to whatever mismatches allowed: is the number oif reads with the iteration mismatch
            
        */
        
        let mut tmp = vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0];



        if ! single_read_input{
            // this is for r1 only if paired end
            for qs in paired_read_quality_score.chars(){
                if qs != '\n'{
                    if qs >= '?'{
                        sample_statistics[sample_id][0] += 1;
                        tmp[0] += 1;
                    }
                    sample_statistics[sample_id][6] += qs as u64 - 33;
                    tmp[6] += qs as u64 - 33;
                    //println!("R1: {}  -> {} -> {}", qs, qs as u64, qs as u64 - 33);
                }
                
            }
            sample_statistics[sample_id][3] += paired_read_quality_score.len() as u64 - 1; 
            tmp[3] += paired_read_quality_score.len() as u64 - 1;   
        }

        tmp_cntr = 0;
        for qs in read_barcode_quality_score.chars(){
            if qs != '\n'{
                if tmp_cntr >= read_barcode_quality_score.len() - barcode_length - 1 {
                    if qs >= '?' {
                        // this is for r3 or barcode
                        sample_statistics[sample_id][2] += 1;
                        tmp[2] += 1;
                    }
                    sample_statistics[sample_id][8] += qs as u64 - 33;
                    tmp[8] += qs as u64 - 33;
                }
                else{
                    // this is for r1 if single end or r2 if paired end
                    if qs >= '?'{
                        sample_statistics[sample_id][shift] += 1;
                        tmp[shift] += 1;
                    }
                    sample_statistics[sample_id][7] += qs as u64 - 33;
                    tmp[7] += qs as u64 - 33;
                
                } 
                //println!("R2/3: {}  -> {} -> {}", qs, qs as u64, qs as u64 - 33);  
                tmp_cntr += 1;
            }
        }
        sample_statistics[sample_id][shift + 3] += (read_barcode_quality_score.len() - barcode_length - 1) as u64; 
        sample_statistics[sample_id][5] += barcode_length as u64;
        
        tmp[shift + 3] += (read_barcode_quality_score.len() - barcode_length - 1) as u64; 
        tmp[5] += barcode_length as u64;
        
        //println!("stat: {} -  {:?}", sample_information[sample_id][0], sample_statistics[sample_id]);
        //println!("This read -  {:?}", tmp);
        
        //if trim_barcode && sample_id < undetermined_label_id{
        
        if read2_has_sequence && trim_barcode{
            read_barcode_seq.truncate(read_barcode_seq.len() - barcode_length - 1);
            read_barcode_seq.push('\n');
            read_barcode_quality_score.truncate(read_barcode_quality_score.len() - barcode_length - 1);
            read_barcode_quality_score.push('\n');
            
        }
        //let dur = start.elapsed();
        //t_trim += dur.as_nanos();


        //let start = Instant::now();
        
        if sample_id == undetermined_label_id
                    && read_cntr > 1000
                    && sample_mismatches[sample_id][0] > (0.75 * read_cntr as f64) as u64
        {
                    panic!( "{}\nAll reads: {}, Undetermined reads: {}",
                         "Seems that there is an issue with the input. Most of the reads are undetermined!",
                         read_cntr, sample_mismatches[undetermined_label_id][0]);
        }
        
        if read2_has_sequence{
            last_piece = out_read_barcode_buffer[sample_id].last_mut().unwrap();
            last_piece.push_str(&read_barcode_header);
            last_piece.push_str(&read_barcode_seq);
            last_piece.push_str(&read_barcode_plus);
            last_piece.push_str(&read_barcode_quality_score);
        }
        
        
        sample_mismatches[sample_id][0] += 1;
        sample_mismatches[sample_id][curr_mismatch + 1] += 1;

        if !single_read_input {
            last_piece = out_paired_read_buffer[sample_id].last_mut().unwrap(); 
            last_piece.push_str(&paired_read_header);
            last_piece.push_str(&paired_read_seq);
            last_piece.push_str(&paired_read_plus);
            last_piece.push_str(&paired_read_quality_score);
        }


        if out_read_barcode_buffer[sample_id].last().unwrap().len() >= read_merging_threshold {
            
            if out_read_barcode_buffer[sample_id].len() >= writing_threshold {
                //let start2 = Instant::now();

                write_reads(&out_read_barcode_buffer[sample_id], &sample_output_path[sample_id].0);
                //let dur2 = start2.elapsed();
                //t_write2 += dur2.as_nanos();
                out_read_barcode_buffer[sample_id].clear();
            }

            out_read_barcode_buffer[sample_id].push(String::new());
        }

        if !single_read_input && out_paired_read_buffer[sample_id].last().unwrap().len() >= read_merging_threshold{
            if out_paired_read_buffer[sample_id].len() >= writing_threshold {
                //let start2 = Instant::now();
                write_reads(&out_paired_read_buffer[sample_id], &sample_output_path[sample_id].1);
                //let dur2 = start2.elapsed();
                //t_write2 += dur2.as_nanos();
                out_paired_read_buffer[sample_id].clear();
            }
            out_paired_read_buffer[sample_id].push(String::new());
                               
        }

        //let dur = start.elapsed();
        //t_write += dur.as_nanos();

        read_cntr += 1;
        
    }
    dur = start.elapsed();
    
    println!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());

    let max_mismatches = allowed_mismatches + 1;
    
    
    for sample_id_itr in 0..out_read_barcode_buffer.len() {
        //println!("{}  ->   {}", sample_id, buffer_info.1.len());
        if out_read_barcode_buffer[sample_id_itr][0].len() > 0 {
            //let start2 = Instant::now();
            write_reads(&out_read_barcode_buffer[sample_id_itr], &sample_output_path[sample_id_itr].0);
            //let dur2 = start2.elapsed();
            //t_write2 += dur2.as_nanos();
            //buffer_info.clear();
        }
        if !single_read_input && out_paired_read_buffer[sample_id_itr][0].len() > 0 {
            write_reads(&out_paired_read_buffer[sample_id_itr], &sample_output_path[sample_id_itr].1);
        }
    }
    //let dur = start.elapsed();
    //t_write += dur.as_nanos();


    //let start = Instant::now();
        
    
    if sample_mismatches[ambiguous_label_id][0] == 0{
        sample_information.pop();
        sample_statistics.pop();
        sample_mismatches.pop();
    }
    


    let mut report_path_main = String::from(&flowcell);
    report_path_main.push('.');
    report_path_main.push_str(&lane);
    report_path_main.push_str(&".mgikit.");
    let mut outfile;
    
    let mut output_file_path;
    let mut report_path;
    
    for sample_index in 0..sample_mismatches.len() {
        sample_statistics[sample_index][9] = sample_mismatches[sample_index][0];
        for cnt in 1..(max_mismatches + 1) {
            sample_statistics[sample_index].push(sample_mismatches[sample_index][cnt]);
        }       
    }
    
    for (project_id, samples) in project_samples.iter() {
        report_path = report_dir_local.clone();
        if project_id != "."{
            println!("generating report for job_number: {} with {} samples.", project_id, samples.len());
            report_path.push_str(&project_id);
            report_path.push('_');
            
        }else{
            println!("generating report for the whole run with {} samples.", sample_mismatches.len());
            
        }
        
        report_path.push_str(&report_path_main);
        
        //println!("{}  -> {:?}", project_id, samples);
        //Start writing info report
        write_index_info_report(&sample_information, &sample_mismatches, samples, max_mismatches, format!("{}{}", &report_path, &"info"));
        //Finish writing info report
        
        //start writing general report
        write_general_info_report(
            &sample_information, &sample_statistics, samples, &flowcell, &lane, format!("{}{}", &report_path, &"general")
        );

    }


    report_path = report_dir_local.clone();
    report_path.push_str(&report_path_main);
    report_path.push_str(&"sample_stats");
    output_file_path = Path::new(&report_path);
    
    let mut out_str = String::from(format!("job_number\tsample_id\tr1_qc_30\tr2_qc_30\tr3_qc_30\tr1_bases\tr2_bases\tr3_bases\tr1_qc\tr2_qc\tr3_qc\tall_reads"));
    for cnt in 0..max_mismatches {
        out_str.push('\t');
        out_str.push_str(&cnt.to_string());
        out_str.push_str(&"-mismatches");
    }
    out_str.push('\n');
    for sample_index in 0..sample_statistics.len() {
        out_str.push_str(&sample_information[sample_index][PROJECT_ID_COLUMN]);
        out_str.push('\t');
        out_str.push_str(&sample_information[sample_index][SAMPLE_COLUMN]);
        for cnt in 0..sample_statistics[0].len(){
            out_str.push('\t');
            out_str.push_str(&sample_statistics[sample_index][cnt].to_string());    
        }
        out_str.push('\n');
    } 
    outfile = File::create(output_file_path).expect("couldn't create output");
    outfile.write_all(&out_str.as_bytes()).unwrap();
    
    
    
    let mut rep_itr = 0;
    report_path = report_dir_local.clone();
    report_path.push_str(&report_path_main);
    report_path.push_str(&"ambiguous_barcode");
    let mut ambiguous_barcodes_out: Vec<_> = ambiguous_barcodes.iter().collect();
    if ambiguous_barcodes_out.len() > 0 {
        output_file_path = Path::new(&report_path);
        outfile = File::create(output_file_path).expect("couldn't create output");
        ambiguous_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
        rep_itr = 0;
        for barcode in &ambiguous_barcodes_out {
            outfile
                .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                .unwrap();
            rep_itr += 1;
            if rep_itr == report_limit {
                break;
            }
        }

        report_path.push_str(".complete");
        output_file_path = Path::new(&report_path);
        outfile = File::create(output_file_path).expect("couldn't create output");
        for barcode in &ambiguous_barcodes_out {
            outfile
                .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                .unwrap();
        }
    }


    //return;
    
    report_path = report_dir_local.clone();
    report_path.push_str(&report_path_main);
    report_path.push_str(&"undetermined_barcode");
    let mut undetermined_barcodes_out: Vec<_> = undetermined_barcodes.iter().collect();
    if undetermined_barcodes.len() > 0 {
        output_file_path = Path::new(&report_path);
        outfile = File::create(output_file_path).expect("couldn't create output");

        undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
        for barcode in &undetermined_barcodes_out {
            outfile
                .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                .unwrap();
            rep_itr += 1;
            if rep_itr == report_limit {
                break;
            }
        }

        report_path.push_str(".complete");
        output_file_path = Path::new(&report_path);
        outfile = File::create(output_file_path).expect("couldn't create output");
        undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
        for barcode in &undetermined_barcodes_out {
            outfile
                .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                .unwrap();
        }
    }

}

/* 
pub fn quality_control_info(read_barcode_file_path: &String, paired_read_file_path: &String){
    //let mut paired_read_data;
    let mut read_barcode_data = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path)).expect("Could not open the file")
    ));
    
    /*if paired_read_file_path.len() > 0 {
        paired_read_data = BufReader::new(MultiGzDecoder::new(
            File::open(Path::new(&paired_read_file_path)).expect("Could not open the file")
        ))
    }*/
    let mut qc_r2_30:u64 = 0;
    let mut qc_r2:f64 = 0.0;
    let mut qc_r1_30:u64 = 0;
    let mut qc_r1:f64 = 0.0;
    let mut r2_bases:u64 = 0;
    let mut curr_line : String;
    let mut read_bytes;
    loop {
        curr_line = String::new();
        
        read_bytes = read_barcode_data.read_line(&mut curr_line).unwrap();
        if read_bytes == 0 {
            //println!("bytes: {}", read_bytes);
            break;
        }

        
        curr_line = String::new();
        read_bytes = read_barcode_data.read_line(&mut curr_line).unwrap();
        
        curr_line = String::new();
        read_bytes = read_barcode_data.read_line(&mut curr_line).unwrap();
        
        curr_line = String::new();
        read_bytes = read_barcode_data.read_line(&mut curr_line).unwrap();
        
        for qs in curr_line.chars(){
            if qs != '\n' {
                if qs >= '?'{
                    qc_r2_30 += 1;
                } 
                qc_r2 += qs as u64 - 33;   
                r2_bases += 1;
            }
        }




    }
    
}
*/
pub fn merge_qc_reports(qc_report_paths: &[String], output_dir: &String){
    if qc_report_paths.len() == 0 {
        panic!("report directories are not provided!");
    }

    if output_dir.len() == 0 {
        println!("output will be written to the current work directory!");
    }

    let mut project_sample_stats_dic : HashMap<String, HashMap<String, Vec<u64>>> = HashMap::new();
    let mut run_dic:HashMap<String, Vec<u64>> = HashMap::new();
    let mut flowcell_id:String = String::new();
    let mut sample_information: Vec<Vec<String>>;
    let mut sample_mismatches:Vec<Vec<u64>>;
    let mut sample_stat_simple:Vec<Vec<u64>>;
    let mut max_mismatches = 0;
    let mut file_content;
    //let mut lines;
    
    let empty_vec = Vec::new();
    let mut path;
    let mut tmp : Vec<&str> ;  
    let mut project_id: String;
    let mut sample_id: String;
    let mut vals  : Vec<String> ;
    
    //let mut lines;
    let mut temp: HashMap<String, Vec<u64>>;
    //let mut sample_list:Vec<String> = Vec::new();
    for qc_report_path in qc_report_paths{
        println!("Reading {} ...", qc_report_path);
        path = Path::new(qc_report_path);
        tmp = path.file_name().unwrap().to_str().unwrap().split(".").collect();
        flowcell_id = tmp[0].to_string();
        file_content = fs::read_to_string(path).unwrap();
        let lines = file_content.lines();
        //let mut header = lines.first().to_string(); 
        //if !lines[0].starts_with("project_id\tsample_id"){
        //    panic!("File does not have th ecorrect format!");
        //}
        for line in lines{
            if line.starts_with("job_number\tsample_id"){
                continue;
            }
            //println!("{}", line);
            vals= line.split("\t").map(|x| x.to_string()).collect();
            project_id = vals.remove(0);
            sample_id = vals.remove(0);
            //sample_list.push(sample_id.to_owned());
            if project_id == "."{
                match run_dic.get_mut(&sample_id){
                    Some(sample_stat) => {
                        for i in 0..sample_stat.len(){
                            sample_stat[i] += vals[i].parse::<u64>().unwrap();
                        }
                    },
                    None => {
                        run_dic.insert(sample_id.to_owned(), vals.into_iter().map(|x| x.parse::<u64>().unwrap()).collect());
                    }
                };
                        
            }else{
                match project_sample_stats_dic.get_mut(&project_id){
                    Some(sample_stats_dic) => {
                        match sample_stats_dic.get_mut(&sample_id){
                            Some(sample_stat) => {
                                for i in 0..sample_stat.len(){
                                    sample_stat[i] += vals[i].parse::<u64>().unwrap();
                                }
                            },
                            None => {
                                sample_stats_dic.insert(sample_id.to_owned(), 
                                vals.into_iter().map(|x| x.parse::<u64>().unwrap()).collect());
                            }
                        };
                    },
                    None => {
                        temp = HashMap::new();
                        temp.insert(sample_id.to_owned(), 
                        vals.into_iter().map(|x| x.parse::<u64>().unwrap()).collect());
                        project_sample_stats_dic.insert(project_id.to_owned(), temp.to_owned());
                    }
                };
            }

            
        }
    }

    
    // copy all samples to the run dictionary
    for (project_id, samples_dic) in project_sample_stats_dic.iter(){
        if project_id == "."{
            panic!("The run should not appear here!");
        }

        for (sample_id, sample_info) in samples_dic.iter(){
            match run_dic.get(sample_id){
                Some(_) => {
                    panic!("Sample ({}) is assigned to a project and the whole run!", sample_id);
                },
                None => {
                    run_dic.insert(sample_id.to_owned(), sample_info.to_owned());
                }
            };
                
        }

    } 

    project_sample_stats_dic.insert(".".to_string(), run_dic);
    let mut sample_info: Vec<(&String, &Vec<u64>)>;
    let mut output_file;
    let lane = String::from("all");
    let mut sample_itr: u64;
    for (project_id, sample_info_dic) in project_sample_stats_dic.iter(){
        
        sample_information = Vec::new();
        sample_mismatches  = Vec::new();
        sample_stat_simple = Vec::new();
    

        sample_info = sample_info_dic.iter().collect();
        //sample_info.sort_by(|a, b| b.1[7].cmp(&a.1[7]));
        
        if project_id == "."{
            output_file = String::new();
        }else{
            output_file = String::from(project_id);
            output_file.push('_');
        }
        output_file.push_str(&flowcell_id);
        output_file.push_str(&".all.mgikit.");
        
        //write_index_info_report(sample_information:&Vec<Vec<String>>, sample_mismatches:&Vec<Vec<u64>>, kept_samples: &Vec<usize>, max_mismatches:usize, output_file:String)
        sample_itr = 0;
        for sample in sample_info{
            
            //println!("{} : {}  -> {:?}", project_id, sample.0, sample.1);
            sample_information.push(vec![sample.0.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), project_id.clone()]);
            sample_mismatches.push(sample.1[7..sample.1.len()].to_owned());
            sample_mismatches.last_mut().unwrap().push(sample_itr);
            sample_stat_simple.push(sample.1[0..9].to_owned());
            max_mismatches = sample.1.len() - 8;
            sample_itr +=1; 
        }
        println!("Writing reports to {} ...", format!("{}/{}(info,general)", output_dir, output_file));
        write_index_info_report(&sample_information, &sample_mismatches, &empty_vec, max_mismatches, 
            format!("{}/{}info", output_dir, output_file));
        write_general_info_report(&sample_information, &sample_stat_simple, &empty_vec, &flowcell_id, &lane, 
            format!("{}/{}general", output_dir, output_file));        
    }


    //for (sample_id, val) in map.iter_mut() {  }
    
}


