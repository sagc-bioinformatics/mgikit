#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::collections::{HashMap,HashSet};
use std::error::Error;
use std::path::Path;
use std::fs::File;
use std::time::Instant;
use std::io::{self, BufRead, Read};
use std::fs::OpenOptions;
use std::io::{BufReader, Write, BufWriter};
use chrono::prelude::*;
use std::fs;
use noodles_bgzf as bgzf;
use memchr::memchr;
use libdeflater::{Compressor, CompressionLvl};
use flate2::read::MultiGzDecoder;


// my modules
mod variables;
use crate::variables::*;
mod samplesheet_utils;
use crate::samplesheet_utils::*;
mod sequence_utils;
use crate::sequence_utils::*;
mod index_dic;
use crate::index_dic::*;



const BUFFER_SIZE: usize = 1 << 19;
//const BLOCK_SIZE: usize = 1 << 16;
//const BUFFER_SIZE_MIN:usize = 1000000;



const HEADER_TAIL: [u8; 5] = [b':', b'N', b':', b'0', b':'];

fn check_file(path_str: &String) {
    if !Path::new(path_str).exists() {
        panic!("File is not accessible: {}", path_str);
    }
}

fn create_folder(path_str: &String) {
    println!("A directory is created: {}", path_str);
    fs::create_dir_all(path_str).unwrap();
}

fn write_illumina_header(output_buffer: &mut [u8], mut buffer_end: usize, mgi_read_header: &[u8], umi: &[u8], l_position:usize) -> usize {
    
    
    output_buffer[buffer_end] = mgi_read_header[l_position + 1];
    buffer_end += 1;

    output_buffer[buffer_end] = b':';
    buffer_end += 1;
    
    for i in l_position + 10..mgi_read_header.len() - 3{
        if  mgi_read_header[i] != b'0'{
            output_buffer[buffer_end..buffer_end + mgi_read_header.len() - 3 - i].copy_from_slice(&mgi_read_header[i..mgi_read_header.len() - 3]);
            buffer_end +=  mgi_read_header.len() - 3 - i;
            break;
        }
    }
     
    output_buffer[buffer_end] = b':';
    buffer_end += 1;

    if mgi_read_header[l_position + 3] != b'0'{
        output_buffer[buffer_end] = mgi_read_header[l_position + 3];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 4];
        output_buffer[buffer_end + 2] = mgi_read_header[l_position + 5];
        buffer_end += 3;
    }else if mgi_read_header[l_position + 4] != b'0'{
        output_buffer[buffer_end] = mgi_read_header[l_position + 4];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 5];
        buffer_end += 2;
    }else{
        output_buffer[buffer_end] = mgi_read_header[l_position + 5];
        buffer_end += 1;
    }
    
    output_buffer[buffer_end] = b':';
    buffer_end += 1;
    
    if mgi_read_header[l_position + 7] != b'0'{
        output_buffer[buffer_end] = mgi_read_header[l_position + 7];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 8];
        output_buffer[buffer_end + 2] = mgi_read_header[l_position + 9];
        buffer_end += 3;
    }else if mgi_read_header[l_position + 8] != b'0'{
        output_buffer[buffer_end] = mgi_read_header[l_position + 8];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 9];
        buffer_end += 2;
    }else{
        output_buffer[buffer_end] = mgi_read_header[l_position + 9];
        buffer_end += 1;
    }


    if umi.len() > 0{
        output_buffer[buffer_end..buffer_end + umi.len()].copy_from_slice(&umi[..]);
        buffer_end += umi.len();
    }
    
    

    output_buffer[buffer_end] = b' ';
    buffer_end += 1;
    
    output_buffer[buffer_end] = mgi_read_header[mgi_read_header.len() - 2];
    buffer_end += 1;
    

    output_buffer[buffer_end..buffer_end + 5].copy_from_slice(&HEADER_TAIL[..]);
    buffer_end += 5;
    buffer_end

}

pub fn write_general_info_report(sample_information:&Vec<Vec<String>>, 
    sample_statistics:&Vec<Vec<u64>>, 
    kept_samples: &Vec<usize>, run:&String, lane:&String, output_file:String){
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
    let non_sample = vec!["undetermined".to_string(), "ambiguous".to_string()];
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
    let mut execluded_reads: f64 = 0.0;
    
    
    let mut unique_sample_ids = HashSet::new();
    for sample_inf in sample_information {
        unique_sample_ids.insert(sample_inf[SAMPLE_COLUMN].clone());
    }

    if unique_sample_ids.len() == sample_information.len(){
        sample_statistics_copy.sort_by(|a, b| b[9].cmp(&a[9]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap())));
    }
    
    
    for sample_id_itr in 0..sample_statistics_copy.len(){
        sample_id = sample_statistics_copy[sample_id_itr].last().unwrap().clone() as usize;
        if kept_samples.len() > 0 && !kept_samples.contains(&sample_id){
            continue;
        }
        
        
        curr_sample_stat = sample_statistics_copy[sample_id_itr].iter().map(|x| x.clone() as f64).collect();

        let mut curr_sample_id: String = sample_information[sample_id][SAMPLE_COLUMN].to_string(); //.split("_").map(|x| x.to_string()).collect::<Vec<String>>()[0].to_lowercase()

        out_str.push_str(&curr_sample_id);
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
        
        curr_sample_id = curr_sample_id.to_lowercase();
        if non_sample.contains(&curr_sample_id){
            execluded_reads += curr_sample_stat[9] as f64;
        }

        if !exclude_non_sample || !non_sample.contains(&curr_sample_id){
            lane_statistics[4] += curr_sample_stat[10];
            //println!("{}  ->  {}   {}    {}    {}", sample_information[sample_id][SAMPLE_COLUMN], lane_statistics[4], curr_sample_stat[10], lane_statistics[1], execluded_reads);
            
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
    final_out_str.push_str(&format!("{:.3?}", lane_statistics[4] / (lane_statistics[1] - execluded_reads) * 100.0));
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
    let mut unique_sample_ids = HashSet::new();
    for sample_inf in sample_information {
        unique_sample_ids.insert(sample_inf[SAMPLE_COLUMN].clone());
    }

    if unique_sample_ids.len() == sample_information.len(){
        sample_mismatches_copy.sort_by(|a, b| b[0].cmp(&a[0]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap())));
    }
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

fn copy_within_a_slice<T: Clone>(v: &mut [T], from: usize, to: usize, len: usize) {
    if from > to {
        let (dst, src) = v.split_at_mut(from);
        dst[to..to + len].clone_from_slice(&src[..len]);
    } else {
        let (src, dst) = v.split_at_mut(to);
        dst[..len].clone_from_slice(&src[from..from + len]);
    }
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
    read1_file_path: &String,
    read2_file_path: &String,
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
    keep_barcode: bool,
    writing_buffer_size: usize,
    comprehensive_scan: bool,
    undetermined_label: &String,
    ambiguous_label: &String,
    force: bool,
    report_limit: usize,
    read1_file_name_suf: &String,
    read2_file_name_suf: &String,
    info_file: &String,
    reporting_level: usize,
    compression_level: u32,
    dynamic_demultiplexing: bool
    
) -> Result<(), Box<dyn Error>> {
    
    // Validate and prepare input data and parameters.

    let mut single_read_input = false;
    let mut info_file_path = if info_file.len() == 0 {
        String::new()
    }else{
        check_file(info_file);
        info_file.to_string()
    };
    let mut paired_read_file_path_final  = read1_file_path.clone();
    let mut read_barcode_file_path_final = read2_file_path.clone();
    
    let mut instrument = String::new();
    let mut run = String::new();
    let mut lane = String::new();
    
    //let quick = true;
    
    if input_folder_path.len() > 0{
        let entries = fs::read_dir(input_folder_path)?
        .map(|res| res.map(|e| e.path()))
        .collect::<Result<Vec<_>, io::Error>>()?;
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
            }else if path.file_name().unwrap().to_str().unwrap().ends_with("BioInfo.csv") && info_file.len() == 0{
                info_file_path = path.to_str().unwrap().to_string();
            }
        }
    }

    if read_barcode_file_path_final.len() == 0{
        //Single End reads
        read_barcode_file_path_final = paired_read_file_path_final;
        paired_read_file_path_final = String::new();
        single_read_input = true;
        println!("Single end read input was detected!");
        println!("Read with Barcode or R1: {}", read_barcode_file_path_final);
    }else{
        println!("Paired ended read input was detected!");
        println!("Paired read or R1: {}", paired_read_file_path_final);
        println!("Read with Barcode or R2: {}", read_barcode_file_path_final);
    }

    if read_barcode_file_path_final.len() == 0 {
        panic!("Input reads are invalid! check the path {}", read_barcode_file_path_final);
    }else if lane.len() == 0{
        let tmp : Vec<&str> = Path::new(&read_barcode_file_path_final).file_name().unwrap().to_str().unwrap().split("_").collect();
        if tmp.len() > 1{
            lane = tmp[1].to_string().trim().to_string();
            if ! lane.starts_with("L0"){
                lane = String::new();
            }
        }
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


    if info_file_path.len() == 0{
        let tmp_path: &Path = Path::new(&read_barcode_file_path_final);
        let tmp_path = tmp_path.with_file_name("BioInfo.csv");
        if tmp_path.exists() {
            info_file_path = tmp_path.to_str().unwrap().to_string();
        }

    }
    println!("info file: {}", info_file_path);

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
    if report_dir_local.chars().last() != Some('/'){
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
  
    
    println!("Output directory: {}", output_directory);
    println!("Reports directory: {}", report_dir_local);
    println!("Instrumnet: {}", instrument);
    println!("Run: {}", run);
    println!("Lane: {}", lane);
    println!("Comprehensive scan mood: {}", comprehensive_scan);
    println!("Compression level: {}. (0 no compression but fast, 9 best compression but slow.)", compression_level);
    println!("Dynamic read determination: {}.", dynamic_demultiplexing);
    

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
    println!("Output buffer size: {}", writing_buffer_size);
    let uncompressed_buffer_size = writing_buffer_size * 2 + 50000;
    let relaxed_writing_buffer_size = ( writing_buffer_size as f32 * 0.95 ) as usize;
    // parse sample/index file and get all mismatches
    
    let tmp_res = parse_sample_index(Path::new(sample_sheet_file_path), &template, i7_rc, i5_rc).unwrap();
    let all_template_data = tmp_res.0;
    let mut sample_information = tmp_res.1;
    let project_samples = tmp_res.2;
    let writing_samples = tmp_res.3;
    let unique_samples_ids = tmp_res.4;
    let barcode_length: usize = all_template_data[0].6[9];
    
    let undetermined_label_id = sample_information.len();
    let ambiguous_label_id = sample_information.len() + 1;
    sample_information.push(vec![undetermined_label.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), String::from(".")]);
    sample_information.push(vec![ambiguous_label.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), String::from(".")]);

    let mut mismatches_dic_i7: Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>  = Vec::new();
    let mut mismatches_dic_i5: Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>  = Vec::new();
    
    for template_details in &all_template_data {
        mismatches_dic_i7.push(get_all_mismatches(&template_details.1, allowed_mismatches));
        mismatches_dic_i5.push(get_all_mismatches(&template_details.2, allowed_mismatches));
    }
    
    
    //  Get reads information 

    let mut whole_read_barcode_len = 0;
    let mut whole_paired_read_len = 0;
    let mut tmp_str = String::new();
    let barcode_read_length: usize;
    let paired_read_length;
    let mut read2_has_sequence = true;
    let only_plus_r1:bool;
    let header_length_r2:usize;
    let mut header_length_r1: usize = 0;

    let bgzf_format_in: bool;

    {
        let mut reader_barcode_read_tmp:bgzf::Reader<File>  = File::open(Path::new(&read_barcode_file_path_final)).map(bgzf::Reader::new).unwrap();
        match reader_barcode_read_tmp.read_line(&mut tmp_str){
            Ok(_)  => {
                bgzf_format_in = true;
            },
            Err(_) => {
                bgzf_format_in = false;
            }
        }
    }
    tmp_str.clear();
    
    let mut reader_barcode_read_tmp = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path_final)).expect("Could not open the file")));
    


    reader_barcode_read_tmp.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    header_length_r2 = tmp_str.len();
    
    let mut l_position = tmp_str.len() - 1;
    for header_chr in tmp_str.chars().rev() {
        if header_chr == 'L' {
            break;
        }

        l_position -= 1;
    }

    if l_position == 0{
            panic!("Can not find the flowcell id in this header {}!", tmp_str);
    }
    let flowcell = tmp_str[1..l_position].to_string();

    println!("Detected flowcell from the header of the first read is {}.", flowcell);
    
    illumina_header_template_p1.push(':');
    illumina_header_template_p1.push_str(&flowcell);
    illumina_header_template_p1.push(':');
    let illumina_header_template_bytes = illumina_header_template_p1.as_bytes();
    
    tmp_str.clear();
    reader_barcode_read_tmp.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    barcode_read_length = tmp_str.chars().count() - 1;
    
    tmp_str.clear();
    reader_barcode_read_tmp.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    
    let only_plus_r2:bool = tmp_str == "+\n";
    if ! only_plus_r2{
        panic!("Expected read format is not satisified. You can try running demultplex-dynamic command.");
    }
    tmp_str.clear();
    reader_barcode_read_tmp.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    
    let mut reader_paired_read = if !single_read_input {
        Some(BufReader::new(MultiGzDecoder::new(
            File::open(Path::new(&paired_read_file_path_final)).expect("Could not open the file")
        )))
    }else{
        None
    };
    
    match reader_paired_read.as_mut() {
        Some(reader_paired_read_buff) => {
            tmp_str.clear();
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            whole_paired_read_len += tmp_str.len();
            header_length_r1 = tmp_str.len();

            tmp_str.clear();
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            paired_read_length = tmp_str.chars().count() - 1;
            whole_paired_read_len += tmp_str.len();
            
            tmp_str.clear();
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            whole_paired_read_len += tmp_str.len();
            only_plus_r1 = tmp_str == "+\n";
            if ! only_plus_r1{
                panic!("Expected read format is not satisified. You can try running demultplex-dynamic command.");
            }

            
            tmp_str.clear();
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            whole_paired_read_len += tmp_str.len();
            
        },
        None => {
            paired_read_length = 0;
        }
    };
    println!("The length of the read with barcode is: {}", barcode_read_length);
    println!("The length of the paired read is: {}", paired_read_length);
    
    
    
    if barcode_length == barcode_read_length{
        read2_has_sequence = false;
        println!("It is assumed that read 2 contains barcode only without read sequence!");    
    }

    
    let barcode_read_actual_length = (barcode_read_length - barcode_length) as u64;
    let paired_read_length_64 = paired_read_length as u64;
    //let barcode_read_length_64 = barcode_read_length as u64;
    let barcode_length_u64 = barcode_length as u64;

    let mut extract_umi;
        
    let mut sample_mismatches: Vec<Vec<u64>> = Vec::new();
    let mut sample_statistics: Vec<Vec<u64>> = Vec::new();
    
    let mut out_read_barcode_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_paired_read_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_read_barcode_buffer_last: Vec<usize> = Vec::new();
    let mut out_paired_read_buffer_last : Vec<usize> = Vec::new();
    
    let mut out_read_barcode_compression_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_paired_read_compression_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_read_barcode_compression_buffer_last: Vec<usize> = Vec::new();
    let mut out_paired_read_compression_buffer_last : Vec<usize> = Vec::new();

    let writen_barcode_length :usize = match trim_barcode {
        true => barcode_length,
        false => 0
    };

    let mut output_barcode_file_writers: Vec<Option<File>> = Vec::new();
    let mut output_paired_file_writers: Vec<Option<File>> = Vec::new();

    for i in 0..sample_information.len(){
        sample_mismatches.push(vec![0; 2 * allowed_mismatches + 2]);
        
        if writing_samples[i] == i {
            let mut tmp = String::from(&output_directory);
            tmp.push_str(&sample_information[i][SAMPLE_COLUMN]);
            
            if i < undetermined_label_id{
                if illumina_format{
                    tmp.push_str(&"_S".to_string());
                    tmp.push_str(&(unique_samples_ids[i] + 1).to_string());
                    tmp.push('_');
                }      
                tmp.push_str(&output_file_format_r2);
            }else{
                if illumina_format{
                    tmp.push('_');
                }
                tmp.push_str(&output_file_format_r2);
            }

            let mut tmp1 = String::from(&output_directory);
            tmp1.push_str(&sample_information[i][SAMPLE_COLUMN]);
            
            if i < undetermined_label_id{
                if illumina_format{
                    tmp1.push_str(&"_S".to_string());
                    tmp1.push_str(&(unique_samples_ids[i] + 1).to_string());
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
            
            if Path::new(&tmp).exists(){
                fs::remove_file(&tmp).unwrap();
            }

            
            let output_file_path = Path::new(&tmp);
            let outfile = OpenOptions::new()
            .append(true)
            .create(true)
            .open(output_file_path)
            .expect("couldn't create output");
            
            output_barcode_file_writers.push(Some(outfile));
            

            if !single_read_input {
                let output_file_path = Path::new(&tmp1);
                let outfile = OpenOptions::new()
                .append(true)
                .create(true)
                .open(output_file_path)
                .expect("couldn't create output");
                output_paired_file_writers.push(Some(outfile));
            }else{
                output_paired_file_writers.push(None);
            }
    
        }else{
            output_barcode_file_writers.push(None);
            output_paired_file_writers.push(None);
        }
        
        
        sample_statistics.push(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        out_read_barcode_buffer_last.push(0);
        out_read_barcode_buffer.push(vec![0; writing_buffer_size]);
        out_read_barcode_compression_buffer_last.push(0);
        out_read_barcode_compression_buffer.push(vec![0; uncompressed_buffer_size]);
        

        out_paired_read_buffer_last.push(0);
        out_paired_read_compression_buffer_last.push(0);
        
        if ! single_read_input{
            out_paired_read_buffer.push( vec![0; writing_buffer_size]);
            out_paired_read_compression_buffer.push( vec![0; uncompressed_buffer_size]);
            
        }

    }

    //let mut curr_template;
    let mut sample_id;
    let mut barcode_read_illumina_header_start:usize = 0;
    let mut barcode_read_illumina_header_end:usize = 0;
    let mut read_cntr:u64 = 0;
    let mut curr_mismatch: usize;
    let mut latest_mismatch:usize;
    let mut curr_umi = String::new();
    let mut curr_barcode = String::new();
    //let mut matching_samples: Vec<(String, u32)> = Vec::new();
    let mut undetermined_barcodes: HashMap<String, u32> = HashMap::new();
    let mut ambiguous_barcodes: HashMap<String, u32> = HashMap::new();
    
    let mut all_i7_mismatches;
    let mut all_i5_mismatches;
    let shift = if single_read_input{0}else{1};

    let start = Instant::now();
    let dur;
    
    

    let mut reader_barcode_read = if bgzf_format_in{
        Box::new(File::open(Path::new(&read_barcode_file_path_final)).map(bgzf::Reader::new).unwrap())
    }else{
        let (tmp_reader, _) = niffler::get_reader(Box::new(std::fs::File::open(read_barcode_file_path_final).unwrap())).unwrap();
        tmp_reader
    };
    
    /* 
    let mut reader_paired_read: Option<Box<dyn Read>> = if !single_read_input {
        if bgzf_format_in{
            Some(File::open(Path::new(&paired_read_file_path_final)).map(bgzf::Reader::new).unwrap())
        }else{
            Some(niffler::get_reader(Box::new(std::fs::File::open(paired_read_file_path_final).unwrap())).unwrap())
        }
    } else {
        None
    };
    
    */
    
    
    //let (mut reader_barcode_read, _)= niffler::get_reader(Box::new(std::fs::File::open(read_barcode_file_path_final).unwrap())).unwrap();
    
    let mut reader_paired_read = if !single_read_input {
        let (s1, _) = niffler::get_reader(Box::new(std::fs::File::open(paired_read_file_path_final).unwrap())).unwrap();
        Some(s1)
    } else {
        None
    };
    

    let mut buffer_1 = [0; BUFFER_SIZE];  // paired read
    let mut buffer_2 = [0; BUFFER_SIZE]; // read with barcode
    let mut curr_bytes:usize;
    
    let mut read_bytes_1: usize = 0;
    if !single_read_input{
        match reader_paired_read{
            Some(ref mut reader) => {
                read_bytes_1 = reader.read(&mut buffer_1[read_bytes_1..]).unwrap();                
            },
            None => panic!("expected single end input!")
        }
    }
    
    let mut read_bytes_2: usize =  reader_barcode_read.read(&mut buffer_2[0..]).unwrap();
            
    if read_bytes_2 == 0{
        panic!("No data!");
    }
        

    let mut header_start: usize = 0;
    let mut read_end: usize;
    
    let mut read_end_pr: usize = 0;
    let mut seq_start: usize;
    let mut plus_start: usize;
    let mut qual_start: usize;
    let mut header_start_pr: usize = 0;
    let mut seq_start_pr: usize = 0;
    let mut plus_start_pr: usize;
    let mut qual_start_pr: usize = 0;
    let mut curr_buffer_end;
    let mut template_itr:usize;
    
    let mut compressor = Compressor::new(CompressionLvl::new(compression_level as i32).unwrap());
    let mut actual_sz:usize;
    let mut curr_writing_sample:usize;
    let total_samples:usize = sample_information.len();
    loop {
        //println!("Read: {}", read_cntr);
        if dynamic_demultiplexing{
            seq_start = match memchr(b'\n', &buffer_2[header_start..read_bytes_2]) {
                None => { panic!("Something wrong with the input data!");},
                Some(loc) => {header_start + loc + 1}
            };
            plus_start = match memchr(b'\n', &buffer_2[seq_start..read_bytes_2]) {
                None => { panic!("Something wrong with the input data!");},
                Some(loc) => {seq_start + loc + 1}
            };
            qual_start = match memchr(b'\n', &buffer_2[plus_start..read_bytes_2]) {
                None => { panic!("Something wrong with the input data!");},
                Some(loc) => {loc + plus_start + 1}
            };
            read_end = match memchr(b'\n', &buffer_2[qual_start..read_bytes_2]) {
                None => { panic!("Something wrong with the input data!");},
                Some(loc) => {qual_start + loc}
            };
        }else{
            seq_start = header_start + header_length_r2;
            plus_start = seq_start + barcode_read_length + 1;
            qual_start = plus_start + 2;
            read_end = barcode_read_length + qual_start;  // \n position.
            if buffer_2[seq_start - 1] != b'\n' || buffer_2[plus_start - 1] != b'\n' ||
                buffer_2[qual_start - 1] != b'\n' || buffer_2[read_end] != b'\n'{
                panic!("Expected read format is not satisified. You can try running demultplex-dynamic command.");
            }     
        }
        
        curr_mismatch = usize::MAX;
        latest_mismatch = usize::MAX;
        sample_id = total_samples;
        template_itr = 0;

        for template_details in &all_template_data {
            let sample_info = &template_details.4;
            let indexes_info = template_details.6;
            
            all_i7_mismatches = &mismatches_dic_i7[template_itr];
            all_i5_mismatches = &mismatches_dic_i5[template_itr];
            
            
            if indexes_info[6] == 1 {
                extract_umi = true;
            } else {
                extract_umi = false;
            }
            
            
            match all_i7_mismatches.get(&buffer_2[(plus_start - indexes_info[2] - 1)
            ..(plus_start - indexes_info[2] + indexes_info[1] - 1)]) {
                Some(i7_matches) => {
                    //println!("{:?}", i7_matches.0);
                    if template_details.5 {
                        match all_i5_mismatches.get(&buffer_2[(plus_start - indexes_info[5] - 1)
                            ..(plus_start - indexes_info[5] + indexes_info[4] - 1)]) {
                            Some(i5_matches) => {
                                //println!("{:?} and {:?}", i7_matches, i5_matches);
                                for i7_match_itr in 0..i7_matches.0.len() {
                                    for i5_match_itr in 0..i5_matches.0.len() {
                                        if latest_mismatch < i7_matches.1 + i5_matches.1{
                                            continue;
                                        }
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
                                                    if sample_id >= total_samples || latest_mismatch > curr_mismatch {
                                                        sample_id = *i5_info;
                                                        latest_mismatch = curr_mismatch;
                                                        curr_barcode = unsafe {
                                                            String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[2] - 1)
                                                    ..(plus_start - indexes_info[2] + indexes_info[1] - 1)].to_vec())
                                                    };
                                                        curr_barcode.push('+');
                                                        curr_barcode.push_str(&unsafe {
                                                            String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[5] - 1)
                                                            ..(plus_start - indexes_info[5] + indexes_info[4] - 1)].to_vec())
                                                        });
            
                                                        if illumina_format {
                                                            if extract_umi {
                                                                curr_umi = String::from(":");
                                                                if i7_rc {
                                                                    curr_umi.push_str(
                                                                        &reverse_complement(
                                                                            &unsafe {
                                                                                String::from_utf8_unchecked(buffer_2[(plus_start  
                                                                                - indexes_info[8]
                                                                                - 1)
                                                                                ..(plus_start
                                                                                    - indexes_info
                                                                                        [8]
                                                                                    + indexes_info
                                                                                        [7]
                                                                                    - 1)].to_vec()
                                                                        )},
                                                                        )
                                                                        .unwrap(),
                                                                    );
                                                                } else {
                                                                    curr_umi.push_str(
                                                                        &unsafe {
                                                                            String::from_utf8_unchecked(buffer_2[(plus_start 
                                                                            - indexes_info[8]
                                                                            - 1)
                                                                            ..(plus_start 
                                                                                - indexes_info[8]
                                                                                + indexes_info[7]
                                                                                - 1)].to_vec())}
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
                                
                            },
                            None => {
                                
                            }
                        }
                    } else {
                        if latest_mismatch < i7_matches.1{
                            continue;
                        }else if latest_mismatch > i7_matches.1{
                            curr_mismatch = i7_matches.1;
                            sample_id = match template_details.4.get(i7_matches.0[0]){
                                Some(tmp) => tmp.0,
                                None => {panic!("There should be a value here");}
                            };
                            //barcode_length = indexes_info[9];
                            if illumina_format {
                                curr_barcode = unsafe {
                                    String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[2] - 1)
                                    ..(plus_start - indexes_info[2] + indexes_info[1] - 1)].to_vec())
                                    };
                                if extract_umi {
                                    curr_umi = String::from(":");
                                    curr_umi.push_str(
                                        &unsafe {
                                            String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[8] - 1)
                                            ..(plus_start - indexes_info[8]
                                                + indexes_info[7]
                                                - 1)].to_vec())}
                                    );
                                    if i7_rc {
                                        curr_umi.push_str(
                                            &reverse_complement(
                                                &unsafe {
                                                    String::from_utf8_unchecked(buffer_2[(plus_start 
                                                    - indexes_info[8]
                                                    - 1)
                                                    ..(plus_start - indexes_info[8]
                                                        + indexes_info[7]
                                                        - 1)].to_vec())}
                                                    
                                            )
                                            .unwrap(),
                                        );
                                    }
                                }
                            }
                        }else{
                            sample_id = ambiguous_label_id;   
                        }
                       
                    }
                },
                None => {
                    
                }
            };
            
            if (sample_id != undetermined_label_id && sample_id < total_samples) && !comprehensive_scan {
                break;
            }
            
            template_itr += 1;
        }
    

        if sample_id >= total_samples {
            sample_id = undetermined_label_id;
            
        }
        
        if sample_id == undetermined_label_id || sample_id == ambiguous_label_id{
            
            if read_cntr > 1000
                    && sample_mismatches[sample_id][0] > (0.75 * read_cntr as f64) as u64
            {
                        panic!( "{}\nAll reads: {}, Undetermined reads: {}",
                                "Seems that there is an issue with the input. Most of the reads are undetermined!",
                                read_cntr, sample_mismatches[undetermined_label_id][0]);
            }
            curr_mismatch = 0;
            if all_template_data.len() == 1
            {
                let indexes_info = all_template_data[0].6;
                curr_barcode = unsafe {
                    String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[2] - 1)
                ..(plus_start - indexes_info[2] + indexes_info[1] - 1)].to_vec())
                };
                
                if  indexes_info[3] == 1{
                    // i5 exist
                    curr_barcode.push('+');
                    curr_barcode.push_str(&unsafe {
                        String::from_utf8_unchecked(buffer_2[(plus_start - indexes_info[5] - 1)
                        ..(plus_start - indexes_info[5] + indexes_info[4] - 1)].to_vec())
                    });
                }
                
            }else{
                curr_barcode = unsafe {
                    String::from_utf8_unchecked(buffer_2[(plus_start - barcode_length - 1)
                ..(plus_start - 1)].to_vec())};
            }
        }

        if reporting_level > 1 {
            if sample_id == undetermined_label_id {
                match undetermined_barcodes.get_mut(&curr_barcode) {
                    Some(barcode_info) => {
                        *barcode_info += 1;
                    }
                    None => {
                        undetermined_barcodes.insert(curr_barcode.clone(), 1);
                    }
                };
            }else if sample_id == ambiguous_label_id {
                match ambiguous_barcodes.get_mut(&curr_barcode) {
                    Some(barcode_info) => {
                        *barcode_info += 1;
                    }
                    None => {
                        ambiguous_barcodes.insert(curr_barcode.clone(), 1);
                    }
                };
                
            }
        }
        

        if !single_read_input {
            if dynamic_demultiplexing{
                seq_start_pr = match memchr(b'\n', &buffer_1[header_start_pr..read_bytes_1]) {
                    None => { panic!("Something wrong with the input data!");},
                    Some(loc) => {header_start_pr + loc + 1}
                };
                plus_start_pr = match memchr(b'\n', &buffer_1[seq_start_pr..read_bytes_1]) {
                    None => { panic!("Something wrong with the input data!");},
                    Some(loc) => {seq_start_pr + loc + 1}
                };
                qual_start_pr = match memchr(b'\n', &buffer_1[plus_start_pr..read_bytes_1]) {
                    None => { panic!("Something wrong with the input data!");},
                    Some(loc) => {loc + plus_start_pr + 1}
                };
                read_end_pr = match memchr(b'\n', &buffer_1[qual_start_pr..read_bytes_1]) {
                    None => { panic!("Something wrong with the input data!");},
                    Some(loc) => {loc + qual_start_pr}
                };
            }else{
                seq_start_pr = header_start_pr + header_length_r1;
                plus_start_pr = seq_start_pr + paired_read_length + 1;
                qual_start_pr = plus_start_pr + 2;
                read_end_pr = paired_read_length + qual_start_pr;
                if buffer_1[seq_start_pr - 1] != b'\n' || buffer_1[plus_start_pr - 1] != b'\n' ||
                    buffer_1[qual_start_pr - 1] != b'\n' || buffer_1[read_end_pr] != b'\n'{
                    panic!("Expected format is not satisified. You can try running demultplex-dynamic command.");
                }     
            }      
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
        
        if reporting_level > 0 {
            if ! single_read_input{
                // this is for r1 only if paired end
                for &qs in buffer_1[qual_start_pr..read_end_pr].iter(){
                        if qs >= 63 {
                            sample_statistics[sample_id][0] += 1;
                            
                        }
                        sample_statistics[sample_id][6] += qs as u64;
                                            
                }
                   
            }

            
            for &qs in buffer_2[read_end - barcode_length..read_end].iter(){
                if qs >= 63 {
                    // this is for r3 or barcode
                    sample_statistics[sample_id][2] += 1;
                    //tmp[2] += 1;
                }
                sample_statistics[sample_id][8] += qs as u64;
                
            }

            for &qs in buffer_2[qual_start..read_end - barcode_length].iter(){
            
                if qs >= 63{
                    sample_statistics[sample_id][shift] += 1;
                    //tmp[shift] += 1;
                }
                sample_statistics[sample_id][7] += qs as u64;
                
            }

        }

        
       
        
                
                
        sample_mismatches[sample_id][0] += 1;
        sample_mismatches[sample_id][curr_mismatch + 1] += 1;
        curr_writing_sample = writing_samples[sample_id];
        // writing preperation
        curr_buffer_end = out_read_barcode_buffer_last[curr_writing_sample];
        // this works for mgi format and unde and ambig and ilumina with a bit of extr
        if curr_buffer_end + read_end - header_start + illumina_header_template_bytes.len() >= writing_buffer_size{
            
            actual_sz = compressor.gzip_compress(&out_read_barcode_buffer[curr_writing_sample][..curr_buffer_end],                   &mut out_read_barcode_compression_buffer[curr_writing_sample][out_read_barcode_compression_buffer_last[curr_writing_sample]..]).unwrap();
            out_read_barcode_compression_buffer_last[curr_writing_sample] += actual_sz;
            curr_buffer_end = 0;

                if out_read_barcode_compression_buffer_last[curr_writing_sample] >= relaxed_writing_buffer_size {
                    match output_barcode_file_writers[curr_writing_sample]{
                        Some(ref mut curr_writer) => {
                            curr_writer.write_all(&out_read_barcode_compression_buffer[curr_writing_sample][..out_read_barcode_compression_buffer_last[curr_writing_sample]]).unwrap();
                            out_read_barcode_compression_buffer_last[curr_writing_sample] = 0;
                        },
                        
                        None => panic!("expeted a writer, but None found!")
                    };
                } 
                         
        }




        if sample_id >= undetermined_label_id{
            out_read_barcode_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end + read_end - header_start + 1].
                    copy_from_slice(&buffer_2[header_start..read_end + 1]);
        
            
            out_read_barcode_buffer_last[curr_writing_sample] = curr_buffer_end + read_end - header_start + 1;
            
        }else if read2_has_sequence
        {
            if illumina_format{
                // Illumina format write the heder and skip the header for mgi.
                barcode_read_illumina_header_start = curr_buffer_end;
                out_read_barcode_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end + illumina_header_template_bytes.len()].
                copy_from_slice(&illumina_header_template_bytes[..]);
                curr_buffer_end += illumina_header_template_bytes.len();
    
                curr_buffer_end = write_illumina_header(&mut out_read_barcode_buffer[curr_writing_sample], 
                    curr_buffer_end, &buffer_2[header_start..seq_start], 
                                        &curr_umi.as_bytes(), l_position);
               
                out_read_barcode_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end + curr_barcode.len()]
                    .copy_from_slice(&curr_barcode.as_bytes());
                curr_buffer_end += curr_barcode.len();
                barcode_read_illumina_header_end = curr_buffer_end;
                header_start = seq_start - 1; 
            }        
            
                        
            out_read_barcode_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end +  plus_start - header_start - writen_barcode_length - 1].
                copy_from_slice(&buffer_2[header_start..plus_start - writen_barcode_length - 1]);
    
            curr_buffer_end += plus_start - header_start - writen_barcode_length - 1;
                    
            out_read_barcode_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end + read_end - plus_start - writen_barcode_length + 1].
                        copy_from_slice(&buffer_2[plus_start - 1..read_end - writen_barcode_length]);
                
            curr_buffer_end += read_end - writen_barcode_length - plus_start + 2; 
            out_read_barcode_buffer[curr_writing_sample][curr_buffer_end - 1] = b'\n';
            out_read_barcode_buffer_last[curr_writing_sample] = curr_buffer_end;
    
        }
        
        if ! single_read_input{
            curr_buffer_end = out_paired_read_buffer_last[curr_writing_sample];
            /* 
            println!("{} - {} - {} - {} - {} --- {} : {}", writing_buffer_size, relaxed_writing_buffer_size, uncompressed_buffer_size, 
                                    curr_buffer_end, 
                                    out_paired_read_compression_buffer_last[curr_writing_sample],
                                compressor.gzip_compress_bound(out_paired_read_compression_buffer_last[curr_writing_sample]),
                                out_paired_read_compression_buffer[curr_writing_sample][out_paired_read_compression_buffer_last[curr_writing_sample]..].len());
            */
            if curr_buffer_end + read_end_pr - header_start_pr + illumina_header_template_bytes.len() >= writing_buffer_size{                
                actual_sz = compressor.gzip_compress(&out_paired_read_buffer[curr_writing_sample][..curr_buffer_end], &mut out_paired_read_compression_buffer[curr_writing_sample][out_paired_read_compression_buffer_last[curr_writing_sample]..]).unwrap();
                out_paired_read_compression_buffer_last[curr_writing_sample] += actual_sz;
                curr_buffer_end = 0;

                if out_paired_read_compression_buffer_last[curr_writing_sample] >= relaxed_writing_buffer_size {
                    match output_paired_file_writers[curr_writing_sample]{
                        Some(ref mut curr_writer) => {
                            curr_writer.write_all(&out_paired_read_compression_buffer[curr_writing_sample][..out_paired_read_compression_buffer_last[curr_writing_sample]]).unwrap();
                            out_paired_read_compression_buffer_last[curr_writing_sample] = 0;
                        },
                        
                        None => panic!("expeted a writer, but None found!")
                    };
                }             
            }

            if illumina_format && sample_id < undetermined_label_id{
                out_paired_read_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end +  barcode_read_illumina_header_end - barcode_read_illumina_header_start].
                copy_from_slice(&out_read_barcode_buffer[curr_writing_sample][barcode_read_illumina_header_start..barcode_read_illumina_header_end]);
                curr_buffer_end += barcode_read_illumina_header_end - barcode_read_illumina_header_start;
                out_paired_read_buffer[curr_writing_sample][curr_buffer_end - curr_barcode.len() - 6] = buffer_1[seq_start_pr - 2];
                header_start_pr = seq_start_pr - 1;
            }

            
            out_paired_read_buffer[curr_writing_sample][curr_buffer_end..curr_buffer_end + read_end_pr - header_start_pr + 1].
                copy_from_slice(&buffer_1[header_start_pr..read_end_pr + 1]);
                
            out_paired_read_buffer_last[curr_writing_sample] = curr_buffer_end + read_end_pr - header_start_pr + 1;
            
            if read_end_pr + whole_paired_read_len >= read_bytes_1 {      
                if read_end_pr < read_bytes_1 - 1 && header_start_pr > 0 {
                    copy_within_a_slice(&mut buffer_1, read_end_pr + 1, 0, read_bytes_1 - read_end_pr - 1);
                    read_bytes_1 -= read_end_pr + 1;
                    header_start_pr = 0;
                }else if read_end_pr == read_bytes_1 - 1{
                    header_start_pr = 0;
                    read_bytes_1 = 0;
                }else{
                    // this for testing, if it doies not appear, we can take it off.
                    panic!("read end should not be greater than the buffer size!");
                }
                
                match reader_paired_read{
                    Some(ref mut reader) => {
                        curr_bytes = reader.read(&mut buffer_1[read_bytes_1..]).unwrap();
                        read_bytes_1 += curr_bytes;
                        
                    },
                    None => panic!("expected sinle end input!")
                };
            }else{
                header_start_pr = read_end_pr + 1;
            }
            
        }

        if read_end + whole_read_barcode_len >= read_bytes_2{
            if read_end < read_bytes_2 - 1 && header_start > 0 {
                copy_within_a_slice(&mut buffer_2, read_end + 1, 0, read_bytes_2 - read_end - 1);
                read_bytes_2 -= read_end + 1;
                header_start = 0;
            }else if read_end == read_bytes_2 - 1{
                header_start = 0;
                read_bytes_2 = 0;
            }else{
                // this for testing, if it doies not appear, we can take it off.
                panic!("read end should not be greater than the buffer size!");
            }
            curr_bytes = reader_barcode_read.read(&mut buffer_2[read_bytes_2..]).unwrap();
            read_bytes_2 += curr_bytes;
            
            if read_bytes_2 == 0{
                break;
            }
            
        }else{
            header_start = read_end + 1;
        }
    
        read_cntr += 1;
        
    }

    let max_mismatches = allowed_mismatches + 1;
    
    
    
    for sample_id in 0..out_read_barcode_buffer.len() {

        sample_statistics[sample_id][3] = paired_read_length_64 * sample_mismatches[sample_id][0];
        sample_statistics[sample_id][6] -= sample_statistics[sample_id][3] * 33;
        sample_statistics[sample_id][shift + 3] = barcode_read_actual_length * sample_mismatches[sample_id][0];
        sample_statistics[sample_id][5] = barcode_length_u64 * sample_mismatches[sample_id][0];
        sample_statistics[sample_id][7] -= sample_statistics[sample_id][shift + 3] * 33;
        sample_statistics[sample_id][8] -= sample_statistics[sample_id][5] * 33;
        curr_writing_sample = writing_samples[sample_id];
        if curr_writing_sample == sample_id 
        {
            if !single_read_input {
                if out_paired_read_buffer_last[curr_writing_sample] > 0{
                    actual_sz = compressor.gzip_compress(&out_paired_read_buffer[curr_writing_sample][..out_paired_read_buffer_last[curr_writing_sample]], &mut out_paired_read_compression_buffer[curr_writing_sample][out_paired_read_compression_buffer_last[curr_writing_sample]..]).unwrap();
                    out_paired_read_compression_buffer_last[curr_writing_sample] += actual_sz;
                
    
                
                    match output_paired_file_writers[curr_writing_sample]{
                            Some(ref mut curr_writer) => {
                                curr_writer.write_all(&out_paired_read_compression_buffer[curr_writing_sample][..out_paired_read_compression_buffer_last[curr_writing_sample]]).unwrap();
                                out_paired_read_compression_buffer_last[curr_writing_sample] = 0;
                                curr_writer.flush().unwrap();
                            },
                            
                            None => panic!("expeted a writer, but None found!")
                    };
                }
        }             
            
    
            if out_read_barcode_buffer_last[curr_writing_sample] > 0{
                actual_sz = compressor.gzip_compress(                    &out_read_barcode_buffer[curr_writing_sample][..out_read_barcode_buffer_last[curr_writing_sample]],             &mut out_read_barcode_compression_buffer[curr_writing_sample][out_read_barcode_compression_buffer_last[curr_writing_sample]..]).unwrap();
                out_read_barcode_compression_buffer_last[curr_writing_sample] += actual_sz;
                
                    match output_barcode_file_writers[curr_writing_sample]{
                        Some(ref mut curr_writer) => {
                            curr_writer.write_all(&out_read_barcode_compression_buffer[curr_writing_sample][..out_read_barcode_compression_buffer_last[curr_writing_sample]]).unwrap();
                            out_read_barcode_compression_buffer_last[curr_writing_sample] = 0;
                            curr_writer.flush().unwrap();
                        },
                        
                        None => panic!("expeted a writer, but None found!")
                    };
                                
            }
            
            
        } 
    }
    

    

    dur = start.elapsed();
    
    println!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());
   
    let start_logs = Instant::now();
    
    
    if sample_mismatches[ambiguous_label_id][0] == 0{
        let sample_id = ambiguous_label_id;
        if sample_mismatches[sample_id][0] == 0{
            let mut tmp = String::from(&output_directory);
            tmp.push_str(&sample_information[sample_id][SAMPLE_COLUMN]);
            let mut tmp1 = String::from(&output_directory);
            tmp1.push_str(&sample_information[sample_id][SAMPLE_COLUMN]);
            
            if illumina_format{
                    tmp.push('_');
                    tmp1.push('_');
            }
            tmp.push_str(&output_file_format_r2);
            tmp1.push_str(&output_file_format_r1);
            
            fs::remove_file(tmp).unwrap();
            if !single_read_input{
                fs::remove_file(tmp1).unwrap();
            }
        }
        sample_information.pop();
        sample_statistics.pop();
        sample_mismatches.pop();

    }
    
    if sample_mismatches[undetermined_label_id][0] == 0{
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
        
        //println!("{}  -> {:.unwrap()}", project_id, samples);
        //Start writing info report
        write_index_info_report(&sample_information, &sample_mismatches, samples, max_mismatches, format!("{}{}", &report_path, &"info"));
        //Finish writing info report
        
        if reporting_level > 0 {
            write_general_info_report(
                &sample_information, &sample_statistics, samples, &flowcell, &lane, format!("{}{}", &report_path, &"general")
            );
        }
        //start writing general report
        

    }

    if reporting_level > 0 {
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
        
    }
    
    if reporting_level > 1 {
        let mut outfile: BufWriter<File>;
        let mut rep_itr = 0;
        report_path = report_dir_local.clone();
        report_path.push_str(&report_path_main);
        report_path.push_str(&"ambiguous_barcode");
        let mut ambiguous_barcodes_out: Vec<_> = ambiguous_barcodes.iter().collect();
        if ambiguous_barcodes_out.len() > 0 {
            output_file_path = Path::new(&report_path);
            outfile = BufWriter::new(File::create(output_file_path).expect("couldn't create output"));
            ambiguous_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
            rep_itr = 0;
            for barcode in &ambiguous_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
                rep_itr += 1;
                if rep_itr == report_limit {
                    break;
                }
            }

            report_path.push_str(".complete");
            output_file_path = Path::new(&report_path);
            outfile = BufWriter::new(File::create(output_file_path).expect("couldn't create output"));
            for barcode in &ambiguous_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
            }
        }


        //return;
        
        report_path = report_dir_local.clone();
        report_path.push_str(&report_path_main);
        report_path.push_str(&"undetermined_barcode");
        let mut undetermined_barcodes_out: Vec<_> = undetermined_barcodes.iter().collect();
        if undetermined_barcodes.len() > 0 {
            output_file_path = Path::new(&report_path);
            outfile = BufWriter::new(File::create(output_file_path).expect("couldn't create output"));

            undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
            for barcode in &undetermined_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
                rep_itr += 1;
                if rep_itr == report_limit {
                    break;
                }
            }

            report_path.push_str(".complete");
            output_file_path = Path::new(&report_path);
            outfile = BufWriter::new(File::create(output_file_path).expect("couldn't create output"));
            undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
            for barcode in &undetermined_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
            }
        }
    }

    let log_dur = start_logs.elapsed();
    
    println!("Writing all logs and reports took {} secs.", log_dur.as_secs());
    Ok(())
}



pub fn detect_template(
    read1_file_path: &String,
    read2_file_path: &String,
    sample_sheet_file_path: &String,
    output_file: &String,
    testing_reads: usize,
    input_barcode_length:usize,
    add_umi:bool,
    use_popular_template: bool,
    max_umi_length: usize
) {
    // Validate input data
    let mut paired_read_file_path_final  = read1_file_path.clone();
    let mut read_barcode_file_path_final = read2_file_path.clone();
    
    let mut single_read_input = false;
    
    if read_barcode_file_path_final.len() == 0{
        //Single End reads
        read_barcode_file_path_final = paired_read_file_path_final;
        paired_read_file_path_final = String::new();
        single_read_input = true;
        println!("Single end read input was detected!");
        if input_barcode_length == 0{
            panic!("Barcode length should be greater than 0 when processing single end reads");
        }
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


    println!("The paired-end read files are assumed to contain the paired reads in the same order!");

    
    let (sample_sheet_header, sample_list, sample_indexes)  = read_sample_sheet_into_dic(Path::new(sample_sheet_file_path)).unwrap();
    
    //println!("{:?}", sample_list);
    //println!("{:?}", sample_sheet_header);
    //println!("{:?}", sample_barcodes);

    
    let mut tmp_str: String = String::new();
    let barcode_read_length: usize;
    let paired_read_length: usize;
    
    let mut reader_barcode_read = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path_final)).expect("Could not open the file")
    ));
    
    reader_barcode_read.read_line(&mut tmp_str).unwrap();
    tmp_str.clear();
    reader_barcode_read.read_line(&mut tmp_str).unwrap();
    barcode_read_length = tmp_str.chars().count() - 1;
    tmp_str.clear();
    let mut reader_paired_read = if !single_read_input {
        Some(BufReader::new(MultiGzDecoder::new(
            File::open(Path::new(&paired_read_file_path_final)).expect("Could not open the file")
        )))
    }else{
        None
    };
    
    match reader_paired_read.as_mut() {
        Some(reader_paired_read_buff) => {
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            tmp_str.clear();
            reader_paired_read_buff.read_line(&mut tmp_str).unwrap();
            paired_read_length = tmp_str.chars().count() - 1;
        },
        None => {
            paired_read_length = 0;
        }
    };
    
    println!("The length of the read with barcode is: {}", barcode_read_length);
    println!("The length of the paired read is: {}", paired_read_length);
    
    if paired_read_length >= barcode_read_length{
        println!("It is assumed that read2 contains barcode only without read sequence!");    
    }

    let barcode_length: usize;
    if input_barcode_length > 0 {
        barcode_length = input_barcode_length as usize;
    }else{
        barcode_length = barcode_read_length - paired_read_length;
        println!("Barcode length is calculated as the difference between R2 length and R1.");
        let max_barcode_length:usize = match sample_indexes.iter().map(|it| it[0].len() + it[2].len()).max(){
                Some(tmp_max) => tmp_max,
                None => panic!("Sample sheet should have samples!") 
            };
            
            if barcode_length > max_barcode_length + max_umi_length{
                panic!("The difference in read length is {}. It is greater than the the length of indexes and possible UMI {}. You need to prvide barcode length for this run!", barcode_length, max_barcode_length + max_umi_length);
            }
        
    }
    
    println!("Barcode length: {}", barcode_length);



    reader_barcode_read = BufReader::new(MultiGzDecoder::new(
        File::open(Path::new(&read_barcode_file_path_final)).expect("Could not open the file")
    ));

    let mut read_barcode_seq;
    let mut read_bytes;

    
    let mut read_cntr:usize = 0;
    let start = Instant::now();
    let mut matches_stat: HashMap<String, Vec<usize>> = HashMap::new();
    
    loop {
        if  read_cntr >= testing_reads {
            //dur = start.elapsed();
            //println!("{} reads  {} seconds", read_cntr, dur.as_secs());
            break;
        }
        read_barcode_seq = String::new();
        read_bytes = reader_barcode_read.read_line(&mut read_barcode_seq).unwrap();
        if read_bytes == 0 {
            break;
        }
        read_barcode_seq = String::new();
        reader_barcode_read.read_line(&mut read_barcode_seq).unwrap();

        //println!("Read seq: {}", read_barcode_seq);
        read_barcode_seq = read_barcode_seq[(read_barcode_seq.len() - barcode_length - 1)..(read_barcode_seq.len() -1)].to_string();
        //println!(" {} Barcode seq: {}", read_cntr, read_barcode_seq);

        for sample_itr in 0..sample_list.len(){
            let sample_i7 = &sample_indexes[sample_itr][0];
            let sample_i7_rc = &sample_indexes[sample_itr][1];
            let sample_i5 = &sample_indexes[sample_itr][2];
            let sample_i5_rc = &sample_indexes[sample_itr][3];

            find_matches(sample_itr, &mut matches_stat, sample_list.len(), 
                &read_barcode_seq, &sample_i7, &sample_i7_rc, &sample_i5, &sample_i5_rc);
            
            
       }
        
            
        reader_barcode_read.read_line(&mut read_barcode_seq).unwrap();
        reader_barcode_read.read_line(&mut read_barcode_seq).unwrap();
        read_cntr += 1;
        
    }

    let mut sample_reads_final: Vec<Vec<(String, usize)>> = Vec::new();                
    
    for sample_itr in 0..sample_list.len(){
        sample_reads_final.push(Vec::new());
        for (template_str, sample_reads) in &matches_stat {
            if sample_reads[sample_itr] > 0{
                sample_reads_final[sample_itr].push((template_str.to_string(), sample_reads[sample_itr]));
            }
        }
        sample_reads_final[sample_itr].sort_by(|a, b| b.1.cmp(&a.1));
    }
    
    let mut popular_template = String::new();
    let mut popular_template_cnt = 0;

    for (template_str, sample_reads) in &matches_stat {
        let cnt = sample_reads.iter().sum();
        if cnt > popular_template_cnt{
            popular_template = template_str.to_string();
            popular_template_cnt = cnt;
        }
    }
    
    if popular_template_cnt == 0{
        panic!("Something wrong, there is no matches with any sample!");
    }

    if use_popular_template{
        let tmp = get_mgikit_template(&popular_template, add_umi, barcode_length, sample_indexes[0][0].len(), sample_indexes[0][2].len());
        println!("The most frequent template is {} - Appeared {} times!", 
        format!("{} i7_rc is {} and i5_rc is {}", tmp.0, tmp.1, tmp.2), 
        popular_template_cnt);
    }
    
    let mut out_str = String::from("sample_id\ti7\ti5\tjob_number\ttemplate\ti7_rc\ti5_rc\n");

    let mut out_str_full = String::from("sample_id\ti7\ti5\tall-matches\n");
    
    for sample_itr in 0..sample_list.len(){
        out_str.push_str(&sample_list[sample_itr][sample_sheet_header[0]]);
        out_str.push('\t');

        out_str.push_str(&sample_list[sample_itr][sample_sheet_header[1]]);
        out_str.push('\t');
        
        if sample_sheet_header[2] == usize::MAX{
            out_str.push('.');
        }else{
            out_str.push_str(&sample_list[sample_itr][sample_sheet_header[2]]);
        }              
        out_str.push('\t');

        if sample_sheet_header[3] == usize::MAX{
            out_str.push('.');
        }else{
            out_str.push_str(&sample_list[sample_itr][sample_sheet_header[3]]);
        }
        
        out_str.push('\t');

        if use_popular_template {
            if sample_reads_final[sample_itr].len() > 0 && popular_template != sample_reads_final[sample_itr][0].0 {
                println!("Using popular template for sample {} instead of its detected template {}!", 
                    sample_list[sample_itr][sample_sheet_header[0]], 
                    sample_reads_final[sample_itr][0].0);
            }else if  sample_reads_final[sample_itr].len() == 0 {
                    println!("Using popular template for sample {} as no matches were found!", 
                    sample_list[sample_itr][sample_sheet_header[0]]);
            }
            
            let template_info = get_mgikit_template(&popular_template, 
                add_umi, 
                barcode_length, 
        sample_list[sample_itr][sample_sheet_header[1]].len(), 
        sample_list[sample_itr][sample_sheet_header[2]].len());
            out_str.push_str(&template_info.0);
            out_str.push('\t');
            out_str.push_str(&template_info.1);
            out_str.push('\t');
            out_str.push_str(&template_info.2);
        }else{

            let template_info = match sample_reads_final[sample_itr].len() > 0{
                true => get_mgikit_template(&sample_reads_final[sample_itr][0].0, 
                add_umi, 
                barcode_length, 
                sample_indexes[sample_itr][0].len(), 
                sample_indexes[sample_itr][2].len()),
                false => (String::from("No matches with this sample"), String::from("."), String::from("."))
                };
            out_str.push_str(&template_info.0);
            out_str.push('\t');
            out_str.push_str(&template_info.1);
            out_str.push('\t');
            out_str.push_str(&template_info.2);
        }


        out_str.push('\n');


        out_str_full.push_str(&sample_list[sample_itr][sample_sheet_header[0]]);
        out_str_full.push('\t');
        
        out_str_full.push_str(&sample_list[sample_itr][sample_sheet_header[1]]);
        out_str_full.push('\t');
        
        if sample_sheet_header[2] == usize::MAX{
            out_str_full.push('.');
        }else{
            out_str_full.push_str(&sample_list[sample_itr][sample_sheet_header[2]]);
        }              
        out_str_full.push('\t');

        if sample_reads_final[sample_itr].len() > 0{
            let tmp: Vec<String> = sample_reads_final[sample_itr].iter().map(|tmp_det| format!("{} ({})", tmp_det.0, tmp_det.1)).collect();
            out_str_full.push_str(&tmp.join("\t"));     
        }else{
            out_str_full.push_str(&String::from("No matches found!"));
        }
        
        out_str_full.push('\n');
    }

    let mut outfile = File::create(Path::new(&format!("{}_template.tsv", output_file))).expect("couldn't create output");
    outfile.write_all(&out_str.as_bytes()).unwrap();

    outfile = File::create(Path::new(&format!("{}_details.tsv", output_file))).expect("couldn't create output");
    outfile.write_all(&out_str_full.as_bytes()).unwrap();

    let dur = start.elapsed();
    
    println!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());

    
}


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
        
        sample_itr = 0;
        for sample in sample_info{
            //println!("{} : {}  -> {:?}", project_id, sample.0, sample.1);
            sample_information.push(vec![sample.0.clone(), String::new(), String::new(), String::new(), String::new(), String::new(), project_id.clone()]);
            sample_mismatches.push(sample.1[9..sample.1.len()].to_owned());
            sample_mismatches.last_mut().unwrap().push(sample_itr);
            sample_stat_simple.push(sample.1[0..11].to_owned());
            max_mismatches = sample.1.len() - 10;
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
