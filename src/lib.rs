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
use memchr::memchr;
use libdeflater::{Compressor, CompressionLvl};
use flate2::read::MultiGzDecoder;
//use std::any::type_name;
//use std::fmt::Debug;
//use glob::glob;

use std::path::PathBuf;

use log::{info, warn};
//use sysinfo::{System, SystemExt};


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
const MAX_SAMPLES: usize = 500;

//const BLOCK_SIZE: usize = 1 << 16;
//const BUFFER_SIZE_MIN:usize = 1000000;



const HEADER_TAIL: [u8; 5] = [b':', b'N', b':', b'0', b':'];

fn check_file<P: AsRef<Path>>(path: &P){
    if !path.as_ref().is_file() {
        panic!("File is not accessible: {}", path.as_ref().display());
    }
}


fn create_folder<P: AsRef<Path>>(path: &P){
    info!("A directory is created: {}", path.as_ref().display());
    fs::create_dir_all(path.as_ref()).unwrap();
}


fn write_illumina_header(output_buffer: &mut [u8], mut buffer_end: usize, mgi_read_header: &[u8], 
    umi: &[u8], l_position:usize, sep_position:usize, sb_header:bool) -> usize {
    output_buffer[buffer_end] = mgi_read_header[l_position + 1];
    buffer_end += 1;

    output_buffer[buffer_end] = b':';
    buffer_end += 1;
    
    
    for i in l_position + 10..sep_position{
        if  mgi_read_header[i] != b'0'{
            output_buffer[buffer_end..buffer_end + sep_position - i].copy_from_slice(&mgi_read_header[i..sep_position]);
            buffer_end +=  sep_position - i;
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

    if sb_header{
        if umi.len() > 0{
            output_buffer[buffer_end] = b':';
            buffer_end += 1;
            output_buffer[buffer_end..buffer_end + umi.len()].copy_from_slice(&umi[..]);
            buffer_end += umi.len();
        }
        output_buffer[buffer_end..buffer_end + mgi_read_header.len() - sep_position - 1]
        .copy_from_slice(&mgi_read_header[sep_position..mgi_read_header.len() - 1]);

        buffer_end += mgi_read_header.len() - sep_position - 1;
        
    }else{
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
    }
    buffer_end

}

fn get_lane_from_file(file_name: &PathBuf, loc: usize, sep: char) -> String {
    let tmp : Vec<&str> = file_name.file_name().unwrap().to_str().unwrap().split(sep).collect();
        if tmp.len() > loc{
            let lane = tmp[loc].to_string();
            if lane.starts_with("L0"){
                return lane;
            }
    }
    String::new()
}

fn validate_and_assigne_input_reads(r1: &String, r2:&String) -> (PathBuf, PathBuf, bool){
    let barcode_read;
    let mut paired_read = String::new();
    
    if r2.len() == 0{
        //Single End reads
        barcode_read = r1.to_string();
        info!("Single end read input was detected!");
        info!("Read with Barcode or R1: {}", barcode_read);
    }else{
        paired_read = r1.clone();
        barcode_read = r2.to_string();
        info!("Paired ended read input was detected!");
        info!("Paired read or R1: {}", paired_read);
        info!("Read with Barcode or R2: {}", barcode_read);
        
        if paired_read.len() == 0 {
            panic!("Input reads are invalid! check the path {}", paired_read);
        }
        check_file(&paired_read);
    }

    if barcode_read.len() == 0 {
        panic!("Input reads are invalid! check the path {}", barcode_read);
    }
    check_file(&barcode_read);
    (PathBuf::from(paired_read), PathBuf::from(barcode_read), r2.len() == 0)
}


fn get_read_files_from_input_dir(in_dir: &String, r1_suf: &String, r2_suf: &String) ->  (PathBuf, PathBuf, bool){
    if in_dir.len() > 0{
        let entries = fs::read_dir(in_dir).expect("can not read directory content!")
        .map(|res| res.map(|e| e.path()))
        .collect::<Result<Vec<_>, io::Error>>().expect("can not collect paths!");
        let mut r1_path: String = String::new();
        let mut r2_path: String = String::new();
        let mut found = 0;
        for path in entries{
            if path.file_name().unwrap().to_str().unwrap().ends_with(r1_suf){
                r1_path = path.to_str().unwrap().to_string();
                found += 1;
            }else if path.file_name().unwrap().to_str().unwrap().ends_with(r2_suf){
                r2_path = path.to_str().unwrap().to_string();
                found += 1;
            }
            if found > 1{
                break;
            }
        }
        if found == 0 {
            panic!("Can not find files that ends with {} or {} under the directory {}", in_dir, r1_suf, r2_suf);
        }
        validate_and_assigne_input_reads(&r1_path, &r2_path)
    }else{
        panic!("input directory is not provided!");
    }
    
}

fn find_info_file<P: AsRef<Path>>(info_file_arg: &P, input_dir: &P, r_file: &P) -> PathBuf{
    if info_file_arg.as_ref().is_file() {
        check_file(info_file_arg);
        return PathBuf::from(info_file_arg.as_ref());
    }else if input_dir.as_ref().is_dir(){
        if let Ok(entries) = fs::read_dir(input_dir) {
            for entry in entries {
                if let Ok(entry) = entry {
                    let file_name = entry.file_name();
                    if file_name.to_string_lossy() == "BioInfo.csv" {
                        return entry.path();
                    }
                }
            }
        }
    }else if r_file.as_ref().is_file() {
        let tmp_path: PathBuf = r_file.as_ref().with_file_name("BioInfo.csv");
        if tmp_path.exists() {
            return tmp_path;
        }
    }
    return PathBuf::new();
}

fn parse_info_file(info_file_path: &PathBuf) -> (String, String){
    let mut instrument: String = String::new();
    let mut run: String = String::new();
    if info_file_path.is_file(){
        let file = File::open(info_file_path).unwrap();
        for line in io::BufReader::new(file).lines() {
            if let Ok(inf) = line {
                //info!("{}", inf);
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
    (instrument, run)
}

fn prepare_output_report_dir(ouput_dir_arg: &String, report_dir_arg: &String, force:bool) -> (PathBuf, PathBuf){
    let use_same_dir :bool;
    let output_directory = if ouput_dir_arg.len() == 0 {
            PathBuf::from(&Local::now().format("mgiKit_%Y%m%dT%H%M%S").to_string())
        }else{
            PathBuf::from(ouput_dir_arg)
        };

    let report_directory = if report_dir_arg.len() == 0 {
            info!("The same output directory will be used for reports.");
            use_same_dir = true;
            output_directory.clone()
        }else{
            use_same_dir = false;
            PathBuf::from(report_dir_arg)
            
        };
    
    
    if output_directory.is_dir(){
        if !force {
            panic!("Output directly exists. Use --froce to overwrite their data: {}", output_directory.display());
        }else{
            info!(
                "Output directory exists. Data will be overwritten at: {}.",
                output_directory.display()
            );
        }
    }else {
        create_folder(&output_directory);
    }

    if report_directory.is_dir(){
        if !use_same_dir{
            if !force {
                panic!("Report directly exists. Use --froce to overwrite their data: {}", report_directory.display());
            }else{
                info!(
                "Report directory exists. Data will be overwritten at: {}.", report_directory.display());
            }
        }
        
    }else {
        create_folder(&report_directory);
    }

    (output_directory, report_directory)
}

fn create_illumina_header_prefix(instrument: &String, run: &String, flowcell:&String) ->String{
    let mut header = String::from("@");
    header.push_str(&instrument);
    header.push(':');
    header.push_str(&run);
    header.push(':');
    header.push_str(&flowcell);
    header.push(':');
    header
}

fn create_output_file_name(sample_name: &String, lane: &String, sample_index: usize, illumina_format:bool) -> (String, String){
    if illumina_format {
        if sample_index == usize::MAX{
            return (format!("{}_{}_R1_001.fastq.gz", sample_name, lane), 
                    format!("{}_{}_R2_001.fastq.gz", sample_name, lane)
                );
        }else{
            return (format!("{}_S{}_{}_R1_001.fastq.gz", sample_name, sample_index, lane), 
                    format!("{}_S{}_{}_R2_001.fastq.gz", sample_name, sample_index, lane)
                );
        }
    }else
    {
        //return (format!("{}_{}_R1.fastq.gz", sample_name, lane), format!("{}_{}_R2.fastq.gz", sample_name, lane));
        return (format!("{}_R1.fastq.gz", sample_name), format!("{}_R2.fastq.gz", sample_name));
    }    
}

fn get_buf_reader(input_file:&PathBuf) -> BufReader<MultiGzDecoder<File>>{
    BufReader::new(MultiGzDecoder::new(
        File::open(input_file).expect("Could not open the file")))
}

fn get_reader(input_file:&PathBuf) -> Box< dyn Read>{
    let (reader, _) = niffler::get_reader(Box::new(std::fs::File::open(input_file).unwrap())).unwrap();
    reader
}

fn get_flowcell_info(mgi_header:&String) ->(String, usize){
    let mut l_position = mgi_header.len() - 1;
    for header_chr in mgi_header.chars().rev() {
        if header_chr == 'L' {
            break;
        }

        l_position -= 1;
    }

    if l_position == 0{
        panic!("Can not find the flowcell id in this header {}!", mgi_header);
    }
    let flowcell = mgi_header[1..l_position].to_string();

    info!("Detected flowcell from the header of the first read is {}.", flowcell);
    (flowcell, l_position)
}

fn parse_sb_file_name(file_name: &String) -> (String, String, String, String){
    //V350170513_L01_23-03891_1.fq.gz
    let parts: Vec<&str> = file_name.split("_").collect();
    (parts[2..parts.len() - 1].join("_"), parts[0].to_string(), parts[1].to_string(), parts[parts.len() - 1].to_string())
}

fn check_output_file(file_path: &PathBuf, force:bool){
    if file_path.exists(){
        if force{
            fs::remove_file(file_path).expect("Can not delete the file");
        }else{
            panic!("Output file exist: {}, you need to use --force to overwrite it", file_path.display());
        }
    }
}

fn hamming_distance(bytes1: &[u8], bytes2: &[u8]) -> Result<usize, &'static str> {
    if bytes1.len() != bytes2.len() {
        return Err("Byte slices must be of equal length");
    }

    Ok(bytes1.iter().zip(bytes2.iter()).filter(|(b1, b2)| b1 != b2).count())
}

pub fn write_general_info_report(sample_information:&Vec<Vec<String>>, 
    sample_statistics:&Vec<Vec<u64>>, 
    kept_samples: &Vec<usize>, run:&String, lane:&String, output_file: &PathBuf){
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

pub fn write_index_info_report(sample_information:&Vec<Vec<String>>, sample_mismatches:&Vec<Vec<u64>>, kept_samples: &Vec<usize>, max_mismatches:usize, output_file: &PathBuf){
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

fn get_read_parts(reader: &mut dyn BufRead) -> (String, String, String, String){
    let mut header: String = String::new();
    let _ = reader.read_line(&mut header).expect("can not read header!");
    let mut seq: String = String::new();
    let _ = reader.read_line(&mut seq).expect("can not read sequence!");
    let mut info: String = String::new();
    let _ = reader.read_line(&mut info).expect("can not read plus!");
    let mut quality: String = String::new();
    let _ = reader.read_line(&mut quality).expect("can not read quality!");

    (header, seq, info, quality)
}

fn read_bytes(reader: &mut Box< dyn Read>, buffer: &mut[u8], minimum:usize, last_byte: &mut usize){
    let mut curr_bytes: usize;
    loop{
        curr_bytes = reader.read(&mut buffer[*last_byte..]).unwrap();
        *last_byte += curr_bytes;
        if *last_byte > minimum || curr_bytes == 0{
            return;
        }
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
    dynamic_demultiplexing: bool,
    compression_buffer_size: usize,
    ignore_undetermined:bool
    
) -> Result<(), Box<dyn Error>> {
    
    // Validate and prepare input data and parameters.
    let trim_barcode = !keep_barcode;
    

    let (paired_read_file_path_final, read_barcode_file_path_final, single_read_input) = 
    if input_folder_path.len() > 0{
        get_read_files_from_input_dir(input_folder_path, read1_file_name_suf, read2_file_name_suf)
    }else{
        validate_and_assigne_input_reads(read1_file_path, read2_file_path)
    };

    let (mut instrument, mut run) = 
                parse_info_file(
            &find_info_file (&PathBuf::from(info_file), 
                                            &PathBuf::new(), 
                                            &read_barcode_file_path_final
                                            )
                                );
    
    
    if arg_instrument.len() > 0{
        instrument = arg_instrument.clone();
    }
    if arg_run.len() > 0{
        run = arg_run.clone();
    }
    
    let lane: String = match arg_lane.len() > 0{
        true => arg_lane.clone(),
        false => get_lane_from_file(&read_barcode_file_path_final , 1, '_')
    };
    
    
    let (output_directory, report_directory) = prepare_output_report_dir(ouput_dir, report_dir, force);
    
    if sample_sheet_file_path.len() == 0 {
        panic!("Sample sheet file is invalid!");
    }
    check_file(sample_sheet_file_path);


    info!("Output directory: {}", output_directory.display());
    info!("Reports directory: {}", report_directory.display());
    info!("Instrumnet: {}", instrument);
    info!("Run: {}", run);
    info!("Lane: {}", lane);
    info!("Comprehensive scan mood: {}", comprehensive_scan);
    info!("Dynamic read determination: {}.", dynamic_demultiplexing);
    info!("Compression level: {}. (0 no compression but fast, 12 best compression but slow.)", compression_level);
    info!("Trim Barcode: {}", trim_barcode);
    // writing_threshold: usize, read_merging_threshold
    info!("Output buffer size: {}", writing_buffer_size);
    if writing_buffer_size < 65536{
        panic!("Writing buffer size '--writing-buffer-size' should be greater than 65536.");
    }
    info!("Compression buffer size: {}", compression_buffer_size);
    if compression_buffer_size > writing_buffer_size{
        panic!("Compression buffer size '--compression-buffer-size' should be less than Writing buffer size ('--writing-buffer-size').");
    }
    info!(
        "Reads that match with multiple samples will be saved in ambiguous_read1/2.fastq file"
    );
    info!("The paired read files are assumed to contain the paired reads in the same order!");
    info!(
        "Allowed mismatches when finding the index are: {}",
        allowed_mismatches
    );

    let trim_barcode = !keep_barcode;
    let illumina_format = !disable_illumina_format;

    
    
    //  Get reads information 

    let mut whole_read_barcode_len: usize ;
    let mut whole_paired_read_len: usize = 0;
    let barcode_read_length: usize;
    let mut paired_read_length: usize = 0;
    let mut read2_has_sequence: bool = true;
    let only_plus_r1:bool;
    let header_length_r2:usize;
    let mut header_length_r1: usize = 0;   
    let flowcell:String;
    let l_position: usize;
    let only_plus_r2: bool;
    
    {
        let mut reader_barcode_read_tmp = get_buf_reader(&read_barcode_file_path_final); 
        let (header, seq, plus, quality) = get_read_parts(&mut reader_barcode_read_tmp);
        
        whole_read_barcode_len = header.len() + seq.len() + plus.len() + quality.len();
        header_length_r2 = header.len();
        (flowcell, l_position) = get_flowcell_info(&header);    
        info!("Detected flowcell from the header of the first read is {}.", flowcell);
        
        barcode_read_length = seq.chars().count() - 1;
        only_plus_r2 = plus == "+\n";
        if ! only_plus_r2 && ! dynamic_demultiplexing{
            panic!("Expected read format is not satisified. You can try rerunning using --flexible parameter.");
        }

        
        if !single_read_input {
            let mut reader_paired_read_buff = get_buf_reader(&paired_read_file_path_final);
            let (header, seq, plus, quality) = get_read_parts(&mut reader_paired_read_buff);
            whole_paired_read_len = header.len() + seq.len() + plus.len() + quality.len();
            header_length_r1 = header.len();
            paired_read_length = seq.chars().count() - 1;
            only_plus_r1 = plus == "+\n";
            if ! only_plus_r1{
                panic!("Expected read format is not satisified. You can try running --flexible command.");
            }
        }
        
    }
    info!("The length of the read with barcode is: {}", barcode_read_length);
    info!("The length of the paired read is: {}", paired_read_length);
    //println!("ZZLENGTH {}  -  {}", whole_paired_read_len, whole_read_barcode_len);
    
    
    let mut illumina_header_prefix_str = String::new();
    let mut illumina_header_prefix = illumina_header_prefix_str.as_bytes();

    if illumina_format {
        info!("Read header and Output files: Illumina format.");
        if lane.len() == 0 || instrument.len() == 0 || run.len() == 0 {
            panic!("Instrument id and run number are required for QC reports and when Illumina format is requested!")
        }
        illumina_header_prefix_str = create_illumina_header_prefix(&instrument, &run, &flowcell);
        illumina_header_prefix = illumina_header_prefix_str.as_bytes(); 
         
    }else {
        info!("Read header and Output files: MGI format.");
        if lane.len() == 0 {
            panic!("Lane number is required for QC reports!")
        }
    }
    
    
    let mut compressor = Compressor::new(CompressionLvl::new(compression_level as i32).unwrap());
    let reqiured_output_buffer_size = writing_buffer_size + compressor.gzip_compress_bound(compression_buffer_size);
    let relaxed_writing_buffer_size = writing_buffer_size; //( writing_buffer_size as f32 * 0.95 ) as usize;
    

    // parse sample/index file and get all mismatches
    let mut i7_rc = arg_i7_rc;
    let i5_rc = arg_i5_rc;
    
    if template.len() > 0 {
        info!(
            "General template is provided and will be used for all samples: {}",
            template
        );
        if i7_rc {
            info!("i7 will be converted to the reverse complement!");
        }

        if i5_rc {
            info!("i5 will be converted to the reverse complement!");
        }
    } else {
        info!("Template will be used from sample/index map file.");
        i7_rc = false;
        //i5_rc = false;
    }

    let tmp_res = parse_sample_index(Path::new(sample_sheet_file_path), &template, i7_rc, i5_rc).unwrap();
    let all_template_data = tmp_res.0;
    let mut sample_information = tmp_res.1;
    let project_samples = tmp_res.2;
    let writing_samples = tmp_res.3;
    let unique_samples_ids = tmp_res.4;
    let barcode_length: usize = all_template_data[0].6[9];
    if barcode_length == barcode_read_length{
        read2_has_sequence = false;
        info!("It is assumed that read 2 contains barcode only without read sequence!");    
    }
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
    
    
    
    let total_samples:usize = sample_information.len();
    let barcode_read_actual_length = (barcode_read_length - barcode_length) as u64;
    let paired_read_length_64 = paired_read_length as u64;
    let barcode_length_u64 = barcode_length as u64;
    let mut extract_umi;

    
    let mut sample_mismatches: Vec<Vec<u64>> = Vec::new();
    let mut sample_statistics: Vec<Vec<u64>> = Vec::new();
    let mut out_read_barcode_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_paired_read_buffer: Vec<Vec<u8>> = Vec::new();
    let mut out_read_barcode_buffer_last: [usize; MAX_SAMPLES] = [0; MAX_SAMPLES];
    let mut out_paired_read_buffer_last : [usize; MAX_SAMPLES] = [0; MAX_SAMPLES];
    
    let writen_barcode_length :usize = match trim_barcode {
        true => barcode_length,
        false => 0
    };

    let mut out_read_barcode_compression_buffer_last: [usize; MAX_SAMPLES] = [0; MAX_SAMPLES];
    let mut out_paired_read_compression_buffer_last : [usize; MAX_SAMPLES] = [0; MAX_SAMPLES];

    let mut extra_comfort_barcode: usize = whole_read_barcode_len + illumina_header_prefix.len() + barcode_length + 25;
    let mut extra_comfort_paired: usize = whole_paired_read_len + illumina_header_prefix.len() + barcode_length + 25;
    

    //let available_memory: u64 = System::new_all().get_available_memory();
    //info!("Available Memory is {} kB", available_memory);

    let mut output_barcode_file_paths: Vec<Option<PathBuf>> = Vec::new();
    let mut output_paired_file_paths: Vec<Option<PathBuf>> = Vec::new();
    let mut output_file_r1: String;
    let mut output_file_r2: String;

    for i in 0..sample_information.len(){
        sample_mismatches.push(vec![0; 2 * allowed_mismatches + 2]);

        if writing_samples[i] == i {
            if i < undetermined_label_id && illumina_format{
                (output_file_r1, output_file_r2) = create_output_file_name(&sample_information[i][SAMPLE_COLUMN], 
                    &lane, unique_samples_ids[i] + 1, true);
            }else if i >= undetermined_label_id && illumina_format{
                (output_file_r1, output_file_r2) = create_output_file_name(&sample_information[i][SAMPLE_COLUMN], 
                    &lane, usize::MAX, true);
            }else{
                (output_file_r1, output_file_r2) = create_output_file_name(&sample_information[i][SAMPLE_COLUMN], 
                    &lane, unique_samples_ids[i] + 1, false);
            }
            
            let barcode_read_output_path = output_directory.join(output_file_r2.clone());
            let paired_read_output_path  = output_directory.join(output_file_r1.clone());

            if barcode_read_output_path.exists(){
                fs::remove_file(&barcode_read_output_path).unwrap();
            }
            
            if ! single_read_input && paired_read_output_path.exists(){
                fs::remove_file(&paired_read_output_path).unwrap();
            }

            
            output_barcode_file_paths.push(Some(barcode_read_output_path));
            output_paired_file_paths.push(Some(paired_read_output_path));
            
            
        }else{
            output_barcode_file_paths.push(None);
            output_paired_file_paths.push(None);
        }
        
        
        sample_statistics.push(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);

        if read2_has_sequence || i >= undetermined_label_id{
            out_read_barcode_buffer.push(vec![0; reqiured_output_buffer_size]);
            
        }else{
            out_read_barcode_buffer.push(vec![]);
        }
        
        out_read_barcode_compression_buffer_last[i] = i * compression_buffer_size;
        if ! single_read_input{
            out_paired_read_buffer.push( vec![0; reqiured_output_buffer_size]);            
            out_paired_read_compression_buffer_last[i] = i * compression_buffer_size;
        }
    }

    let mut out_read_barcode_compression_buffer: Vec<u8> = vec![0; compression_buffer_size * total_samples];

    let mut out_paired_read_compression_buffer = if ! single_read_input{
        vec![0; compression_buffer_size * total_samples]
    }else{
        Vec::new()
    };

    if all_template_data.len() > 1 {
        info!("Mixed library is detected! different barcode templates for some samples!");
    }else{
        info!("Same barcode template is used for all samples!");
    }
    for template_details in &all_template_data {
        let indexes_info = template_details.6;
    }


    
    //let mut curr_template;
    let mut sample_id;
    let mut barcode_read_illumina_header_start:usize = 0;
    let mut barcode_read_illumina_header_end:usize = 0;
    let mut read_cntr:u64 = 1;
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
    
    

    let mut reader_barcode_read = get_reader(&read_barcode_file_path_final); 
    let mut reader_paired_read = if !single_read_input {
        let r1 = get_reader(&paired_read_file_path_final);
        Some(r1)
    } else {
        None
    };
    

    let mut buffer_1 = [0; BUFFER_SIZE];  // paired read
    let mut buffer_2 = [0; BUFFER_SIZE]; // read with barcode
    
    let mut read_bytes_1: usize = 0;
    if !single_read_input{
        match reader_paired_read{
            Some(ref mut reader) => {
                read_bytes(reader, &mut buffer_1, whole_paired_read_len, &mut read_bytes_1);
            },
            None => panic!("expected single end input!")
        }
        if read_bytes_1 == 0{
            panic!("No data in R1!");
        }
    }
        
    let mut read_bytes_2: usize = 0;
    read_bytes(&mut reader_barcode_read, &mut buffer_2, whole_read_barcode_len, &mut read_bytes_2);
    if read_bytes_2 == 0{
        panic!("No data in R2!");
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
    
    let mut actual_sz:usize;
    let mut curr_writing_sample:usize;
    
    let mut curr_writer: File;
    let mut curr_buffer_limit:usize;
    let mut undertmined_threshold_check = 1000;
    loop {
        //info!("Read: {}", read_cntr);
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
            whole_read_barcode_len = 2 * (read_end - header_start + 1);
            extra_comfort_barcode  = whole_read_barcode_len + illumina_header_prefix.len() + barcode_length + 25;
    
            

        }else{
            seq_start = header_start + header_length_r2;
            plus_start = seq_start + barcode_read_length + 1;
            qual_start = plus_start + 2;
            read_end = barcode_read_length + qual_start;  // \n position.
            if buffer_2[seq_start - 1] != b'\n' || buffer_2[plus_start - 1] != b'\n' ||
                buffer_2[qual_start - 1] != b'\n' || buffer_2[read_end] != b'\n'{
                panic!("Expected read format is not satisified. You can try rerunning using --flexible parameter.");
            }
            if read_end >= read_bytes_2{
                panic!("Something wrong, read 2 end is beyond the actual buffer size!.");
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
            
            //println!("{:?}", indexes_info);
            //println!("{}", template_details.5);
            
            if indexes_info[6] == 1 {
                extract_umi = true;
            } else {
                extract_umi = false;
            }
            
            
            match all_i7_mismatches.get(&buffer_2[(plus_start - indexes_info[2] - 1)
            ..(plus_start - indexes_info[2] + indexes_info[1] - 1)]) {
                Some(i7_matches) => {
                    //info!("{:?}", i7_matches.0);
                    if template_details.5 {
                        match all_i5_mismatches.get(&buffer_2[(plus_start - indexes_info[5] - 1)
                            ..(plus_start - indexes_info[5] + indexes_info[4] - 1)]) {
                            Some(i5_matches) => {
                                //info!("{:?} and {:?}", i7_matches, i5_matches);
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
            
            if read_cntr > undertmined_threshold_check
                    && sample_mismatches[sample_id][0] > (0.75 * read_cntr as f64) as u64
            {
                        if ignore_undetermined{
                            warn!( "{}\nAll reads: {}, Undetermined reads: {}",
                            "Seems that there is an issue with the input. Most of the reads are undetermined!",
                            read_cntr, sample_mismatches[undetermined_label_id][0]);
                            undertmined_threshold_check = u64::MAX;
                        }else{
                            panic!( "{}\nAll reads: {}, Undetermined reads: {}",
                                "Seems that there is an issue with the input. Most of the reads are undetermined!",
                                read_cntr, sample_mismatches[undetermined_label_id][0]);
                        }
                        
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
                //println!("Z1: {} - {} - {} - {} - {}", header_start_pr, seq_start_pr, plus_start_pr, qual_start_pr, read_bytes_1);
                read_end_pr = match memchr(b'\n', &buffer_1[qual_start_pr..read_bytes_1]) {
                    None => {

                        println!("{}", unsafe { std::str::from_utf8(&buffer_1[header_start_pr..read_bytes_1]).unwrap()}); 
                        panic!("Something wrong with the input data!");
                    },
                    Some(loc) => {loc + qual_start_pr}
                };
                
                whole_paired_read_len = 2* (read_end_pr - header_start_pr + 1);
                extra_comfort_paired = whole_paired_read_len + illumina_header_prefix.len() + barcode_length + 25;
                //println!("Z2: {} - {}  - {}  -> {}", header_start_pr, read_end_pr, read_cntr, whole_paired_read_len);
                //println!("{}", unsafe { std::str::from_utf8(&buffer_1[header_start_pr..read_end_pr + 1]).unwrap()}); 

            }else{
                seq_start_pr = header_start_pr + header_length_r1;
                plus_start_pr = seq_start_pr + paired_read_length + 1;
                qual_start_pr = plus_start_pr + 2;
                read_end_pr = paired_read_length + qual_start_pr;
                if buffer_1[seq_start_pr - 1] != b'\n' || buffer_1[plus_start_pr - 1] != b'\n' ||
                    buffer_1[qual_start_pr - 1] != b'\n' || buffer_1[read_end_pr] != b'\n'{
                    panic!("Expected format is not satisified. You can try rerunning using --flexible parameter.");
                }
                if read_end_pr >= read_bytes_1{
                    panic!("Something wrong, read 2 end is beyond the actual buffer size!.");
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
        curr_buffer_end = out_read_barcode_compression_buffer_last[curr_writing_sample];
        curr_buffer_limit =  curr_writing_sample * compression_buffer_size + compression_buffer_size;
        // this works for mgi format and unde and ambig and ilumina with a bit of extr
        if curr_buffer_end + extra_comfort_barcode >= curr_buffer_limit {
            actual_sz = compressor.gzip_compress(&out_read_barcode_compression_buffer[curr_writing_sample * compression_buffer_size..curr_buffer_end],
                                   &mut out_read_barcode_buffer[curr_writing_sample][out_read_barcode_buffer_last[curr_writing_sample]..]).unwrap();
            out_read_barcode_buffer_last[curr_writing_sample] += actual_sz;

            curr_buffer_end = curr_writing_sample * compression_buffer_size;

            if out_read_barcode_buffer_last[curr_writing_sample] >= relaxed_writing_buffer_size {
                match &output_barcode_file_paths[curr_writing_sample]{
                    Some(output_file_path) => {
                        curr_writer = OpenOptions::new()
                        .append(true)
                        .create(true)
                        .open(output_file_path)
                        .expect("couldn't create output");
                        curr_writer.write_all(&out_read_barcode_buffer[curr_writing_sample][..out_read_barcode_buffer_last[curr_writing_sample]]).unwrap();
                        curr_writer.flush().expect("couldn't flush output");
                        out_read_barcode_buffer_last[curr_writing_sample] = 0;
                    },
                    
                    None => panic!("expeted a writer, but None found!")
                };
            }          
        }




        if sample_id >= undetermined_label_id{
            out_read_barcode_compression_buffer[curr_buffer_end..curr_buffer_end + read_end - header_start + 1].
                    copy_from_slice(&buffer_2[header_start..read_end + 1]);
        
            
            out_read_barcode_compression_buffer_last[curr_writing_sample] = curr_buffer_end + read_end - header_start + 1;
            
        }else if read2_has_sequence
        {
            if illumina_format{
                // Illumina format write the heder and skip the header for mgi.
                barcode_read_illumina_header_start = curr_buffer_end;
                out_read_barcode_compression_buffer[curr_buffer_end..curr_buffer_end + illumina_header_prefix.len()].
                copy_from_slice(&illumina_header_prefix[..]);
                curr_buffer_end += illumina_header_prefix.len();
    
                curr_buffer_end = write_illumina_header(&mut out_read_barcode_compression_buffer, 
                    curr_buffer_end, &buffer_2[header_start..seq_start], 
                                        &curr_umi.as_bytes(), l_position, seq_start - header_start - 3, false);
               
                out_read_barcode_compression_buffer[curr_buffer_end..curr_buffer_end + curr_barcode.len()]
                    .copy_from_slice(&curr_barcode.as_bytes());
                curr_buffer_end += curr_barcode.len();
                barcode_read_illumina_header_end = curr_buffer_end;
                header_start = seq_start - 1; 
            }        
            
                        
            out_read_barcode_compression_buffer[curr_buffer_end..curr_buffer_end +  plus_start - header_start - writen_barcode_length - 1].
                copy_from_slice(&buffer_2[header_start..plus_start - writen_barcode_length - 1]);
    
            curr_buffer_end += plus_start - header_start - writen_barcode_length - 1;
                    
            out_read_barcode_compression_buffer[curr_buffer_end..curr_buffer_end + read_end - plus_start - writen_barcode_length + 1].
                        copy_from_slice(&buffer_2[plus_start - 1..read_end - writen_barcode_length]);
                
            curr_buffer_end += read_end - writen_barcode_length - plus_start + 2; 
            out_read_barcode_compression_buffer[curr_buffer_end - 1] = b'\n';
            out_read_barcode_compression_buffer_last[curr_writing_sample] = curr_buffer_end;

            if curr_buffer_end >= curr_buffer_limit{
                // remove if it does not happen again!
                panic!("R2 exceeds the buffer limit! somthing wrong with writing to the buffers!");
            }
    
        }
        
        if ! single_read_input{
            curr_buffer_end = out_paired_read_compression_buffer_last[curr_writing_sample];
            if curr_buffer_end + extra_comfort_paired >= curr_buffer_limit{                
                
                actual_sz = compressor.gzip_compress(&out_paired_read_compression_buffer[curr_writing_sample * compression_buffer_size .. curr_buffer_end], 
                    &mut out_paired_read_buffer[curr_writing_sample][out_paired_read_buffer_last[curr_writing_sample]..]).unwrap();
                out_paired_read_buffer_last[curr_writing_sample] += actual_sz;
                curr_buffer_end = curr_writing_sample * compression_buffer_size;

                if out_paired_read_buffer_last[curr_writing_sample] >= relaxed_writing_buffer_size {
                    match &output_paired_file_paths[curr_writing_sample]{
                        Some(output_file_path) => {
                            curr_writer = OpenOptions::new()
                            .append(true)
                            .create(true)
                            .open(output_file_path)
                            .expect("couldn't create output");
                            curr_writer.write_all(&out_paired_read_buffer[curr_writing_sample][..out_paired_read_buffer_last[curr_writing_sample]]).unwrap();
                            curr_writer.flush().expect("couldn't flush output");
                            out_paired_read_buffer_last[curr_writing_sample] = 0;
                        },
                        
                        None => panic!("expeted a writer, but None found!")
                    };
                }             
            }

            if illumina_format && sample_id < undetermined_label_id{
                
                if read2_has_sequence{
                    out_paired_read_compression_buffer[curr_buffer_end..curr_buffer_end +  barcode_read_illumina_header_end - barcode_read_illumina_header_start].
                    copy_from_slice(&out_read_barcode_compression_buffer[barcode_read_illumina_header_start..barcode_read_illumina_header_end]);
                    curr_buffer_end += barcode_read_illumina_header_end - barcode_read_illumina_header_start;
                    out_paired_read_compression_buffer[curr_buffer_end - curr_barcode.len() - 6] = buffer_1[seq_start_pr - 2];
                }else{

                    out_paired_read_compression_buffer[curr_buffer_end..curr_buffer_end + illumina_header_prefix.len()].
                    copy_from_slice(&illumina_header_prefix[..]);
                    curr_buffer_end += illumina_header_prefix.len();
    
                    curr_buffer_end = write_illumina_header(&mut out_paired_read_compression_buffer, 
                        curr_buffer_end, &buffer_1[header_start_pr..seq_start_pr], 
                                            &curr_umi.as_bytes(), l_position, seq_start_pr - header_start_pr - 3, false);
                   
                    out_paired_read_compression_buffer[curr_buffer_end..curr_buffer_end + curr_barcode.len()]
                        .copy_from_slice(&curr_barcode.as_bytes());
                    curr_buffer_end += curr_barcode.len();

                }
                
                header_start_pr = seq_start_pr - 1;
            }

            
            out_paired_read_compression_buffer[curr_buffer_end..curr_buffer_end + read_end_pr - header_start_pr + 1].
                copy_from_slice(&buffer_1[header_start_pr..read_end_pr + 1]);
                
            out_paired_read_compression_buffer_last[curr_writing_sample] = curr_buffer_end + read_end_pr - header_start_pr + 1;
            
            if out_paired_read_compression_buffer_last[curr_writing_sample] >= curr_buffer_limit{
                // remove if it does not happen again!
                panic!("R1 exceeds the buffer limit! somthing wrong with writing to the buffers!");
            }

            
            if read_end_pr + whole_paired_read_len >= read_bytes_1 {
                //println!("R1: {}-{}-{}-{}-{}", read_cntr, whole_paired_read_len, header_start_pr, read_end_pr, read_bytes_1);
                if read_end_pr < read_bytes_1 - 1 {
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
                        read_bytes(reader, &mut buffer_1, whole_paired_read_len, &mut read_bytes_1);
                    },
                    None => panic!("expected sinle end input!")
                };
            }else{
                header_start_pr = read_end_pr + 1;
            }
            
        }

        if read_end + whole_read_barcode_len >= read_bytes_2{
            //println!("R2: {}-{}-{}-{}-{}", read_cntr, whole_read_barcode_len, header_start, read_end, read_bytes_2);
            if read_end < read_bytes_2 - 1 {
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
            
            read_bytes(&mut reader_barcode_read, &mut buffer_2, whole_read_barcode_len, &mut read_bytes_2);
            if read_bytes_2 == 0{
                if !single_read_input && read_bytes_1 != 0 {
                    panic!("Seems like R1 buffer still have reads! It is expected to hav same read count as R2!")
                }
                break;
            }
            
        }else{
            header_start = read_end + 1;
        }
    
        read_cntr += 1;
        
    }

    let max_mismatches = allowed_mismatches + 1;
    
    
    
    for sample_id in 0..total_samples {

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
                if out_paired_read_compression_buffer_last[curr_writing_sample] > curr_writing_sample * compression_buffer_size {
                    actual_sz = compressor.gzip_compress(&out_paired_read_compression_buffer[curr_writing_sample * compression_buffer_size..out_paired_read_compression_buffer_last[curr_writing_sample]],
                         &mut out_paired_read_buffer[curr_writing_sample][out_paired_read_buffer_last[curr_writing_sample]..]).unwrap();
                    out_paired_read_buffer_last[curr_writing_sample] += actual_sz;
                }

                if out_paired_read_buffer_last[curr_writing_sample] > 0{
                    match &output_paired_file_paths[curr_writing_sample]{
                        Some(output_file_path) => {
                            curr_writer = OpenOptions::new()
                            .append(true)
                            .create(true)
                            .open(output_file_path)
                            .expect("couldn't create output");
                            curr_writer.write_all(&out_paired_read_buffer[curr_writing_sample][..out_paired_read_buffer_last[curr_writing_sample]]).unwrap();
                            out_paired_read_buffer_last[curr_writing_sample] = 0;
                            curr_writer.flush().unwrap();
                        },
                        None => panic!("expeted a writer, but None found!")
                    };
                }
            }             
            
    
            if out_read_barcode_compression_buffer_last[curr_writing_sample] > curr_writing_sample * compression_buffer_size{
                actual_sz = compressor.gzip_compress(&out_read_barcode_compression_buffer[curr_writing_sample * compression_buffer_size ..out_read_barcode_compression_buffer_last[curr_writing_sample]],             
                    &mut out_read_barcode_buffer[curr_writing_sample][out_read_barcode_buffer_last[curr_writing_sample]..]).unwrap();
                out_read_barcode_buffer_last[curr_writing_sample] += actual_sz;
            }
            if out_read_barcode_buffer_last[curr_writing_sample] > 0{
                match &output_barcode_file_paths[curr_writing_sample]{
                    Some(output_file_path) => {
                        curr_writer = OpenOptions::new()
                        .append(true)
                        .create(true)
                        .open(output_file_path)
                        .expect("couldn't create output");
                        curr_writer.write_all(&out_read_barcode_buffer[curr_writing_sample][..out_read_barcode_buffer_last[curr_writing_sample]]).unwrap();
                        out_read_barcode_buffer_last[curr_writing_sample] = 0;
                        curr_writer.flush().unwrap();
                    },
                    None => panic!("expeted a writer, but None found!")
                };            
            }
            
            
        } 
    }
    

    

    dur = start.elapsed();
    
    info!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());
   
    let start_logs = Instant::now();
    
    
    if sample_mismatches[ambiguous_label_id][0] == 0{
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
    
    let mut file_name_extra: String ;
    for sample_index in 0..sample_mismatches.len() {
        sample_statistics[sample_index][9] = sample_mismatches[sample_index][0];
        for cnt in 1..(max_mismatches + 1) {
            sample_statistics[sample_index].push(sample_mismatches[sample_index][cnt]);
        }       
    }
    
    for (project_id, samples) in project_samples.iter() {
        
        if project_id != "."{
            info!("generating report for job_number: {} with {} samples.", project_id, samples.len());
            file_name_extra = project_id.clone();
            file_name_extra.push('_');
            
        }else{
            info!("generating report for the whole run with {} samples.", sample_mismatches.len());
            file_name_extra = String::new();
        }
        
        
        //info!("{}  -> {:.unwrap()}", project_id, samples);
        //Start writing info report
        write_index_info_report(&sample_information, &sample_mismatches, samples, max_mismatches, &report_directory.clone().join(format!("{}{}info", &file_name_extra, &report_path_main)));
        //Finish writing info report
        
        if reporting_level > 0 {
            write_general_info_report(
                &sample_information, &sample_statistics, samples, &flowcell, &lane, &report_directory.clone().join(format!("{}{}general", &file_name_extra, &report_path_main)));
        }
        //start writing general report
        

    }

    if reporting_level > 0 {
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
        outfile = File::create(report_directory.clone().join(format!("{}sample_stats", &report_path_main))).expect("couldn't create output");
        outfile.write_all(&out_str.as_bytes()).unwrap();
        
    }
    
    if reporting_level > 1 {
        let mut outfile: BufWriter<File>;
        let mut rep_itr = 0;
        let mut ambiguous_barcodes_out: Vec<_> = ambiguous_barcodes.iter().collect();
        if ambiguous_barcodes_out.len() > 0 {
            outfile = BufWriter::new(File::create(report_directory.clone().join(format!("{}ambiguous_barcode", &report_path_main))).expect("couldn't create output"));
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
            
            outfile = BufWriter::new(File::create(report_directory.clone().join(format!("{}ambiguous_barcode.complete", &report_path_main))).expect("couldn't create output"));
            for barcode in &ambiguous_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
            }
        }


        //return;
        
        let mut undetermined_barcodes_out: Vec<_> = undetermined_barcodes.iter().collect();
        if undetermined_barcodes.len() > 0 {
            outfile = BufWriter::new(File::create(report_directory.clone().join(&format!("{}undetermined_barcode", report_path_main))).expect("couldn't create output"));
            undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
            for barcode in &undetermined_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
                rep_itr += 1;
                if rep_itr == report_limit {
                    break;
                }
            }

            outfile = BufWriter::new(File::create(report_directory.clone().join(&format!("{}undetermined_barcode.complete", report_path_main))).expect("couldn't create output"));
            undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
            for barcode in &undetermined_barcodes_out {
                outfile
                    .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes()).unwrap();
            }
        }
    }

    let log_dur = start_logs.elapsed();
    
    info!("Writing all logs and reports took {} secs.", log_dur.as_secs());
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
        info!("Single end read input was detected!");
        if input_barcode_length == 0{
            panic!("Barcode length should be greater than 0 when processing single end reads");
        }
    }else{
        info!("Paired ended read input was detected!");
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


    info!("The paired-end read files are assumed to contain the paired reads in the same order!");

    
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
    
    info!("The length of the read with barcode is: {}", barcode_read_length);
    info!("The length of the paired read is: {}", paired_read_length);
    
    if paired_read_length >= barcode_read_length{
        info!("It is assumed that read2 contains barcode only without read sequence!");    
    }

    let barcode_length: usize;
    if input_barcode_length > 0 {
        barcode_length = input_barcode_length as usize;
    }else{
        barcode_length = barcode_read_length - paired_read_length;
        info!("Barcode length is calculated as the difference between R2 length and R1.");
        let max_barcode_length:usize = match sample_indexes.iter().map(|it| it[0].len() + it[2].len()).max(){
                Some(tmp_max) => tmp_max,
                None => panic!("Sample sheet should have samples!") 
            };
            
            if barcode_length > max_barcode_length + max_umi_length{
                panic!("The difference in read length is {}. It is greater than the the length of indexes and possible UMI {}. You need to prvide barcode length for this run!", barcode_length, max_barcode_length + max_umi_length);
            }
        
    }
    
    info!("Barcode length: {}", barcode_length);



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
        info!("The most frequent template is {} - Appeared {} times!", 
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
                info!("Using popular template for sample {} instead of its detected template {}!", 
                    sample_list[sample_itr][sample_sheet_header[0]], 
                    sample_reads_final[sample_itr][0].0);
            }else if  sample_reads_final[sample_itr].len() == 0 {
                    info!("Using popular template for sample {} as no matches were found!", 
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
    
    info!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());

    
}


pub fn merge_qc_reports(qc_report_paths: &[String], 
    output_dir: &String,
    lane: &String,
    project: &String){
    if qc_report_paths.len() == 0 {
        panic!("report directories are not provided!");
    }

    if output_dir.len() == 0 {
        info!("output will be written to the current work directory!");
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
        info!("Reading {} ...", qc_report_path);
        path = Path::new(qc_report_path);
        tmp = path.file_name().unwrap().to_str().unwrap().split(".").collect();
        flowcell_id = if project.len() == 0 {
            tmp[0].to_string()
        }else{
            project.to_string()
        };

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
    let used_lane = if lane.len() == 0 {
            String::from("all")
        }else{
            lane.to_string()
        };
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
        output_file.push('.');
        output_file.push_str(&used_lane);
        output_file.push_str(&".mgikit.");
        
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
        info!("Writing reports to {} ...", format!("{}/{}(info,general)", output_dir, output_file));
        write_index_info_report(&sample_information, &sample_mismatches, &empty_vec, max_mismatches, 
            &PathBuf::from(format!("{}/{}info", output_dir, output_file)));
        write_general_info_report(&sample_information, &sample_stat_simple, &empty_vec, &flowcell_id, &used_lane, 
            &PathBuf::from(format!("{}/{}general", output_dir, output_file)));        
    }


    //for (sample_id, val) in map.iter_mut() {  }
    
}


pub fn reformat(
    r1: &String,
    r2: &String,
    ouput_dir: &String,
    report_dir: &String,
    arg_lane: &String,
    arg_instrument: &String,
    arg_run: &String,
    illumina_format:bool,
    writing_buffer_size: usize,
    force: bool,
    report_limit: usize,
    info_file: &String,
    reporting_level: usize,
    compression_level: u32,
    compression_buffer_size: usize,
    umi_length:usize,
    sample_index:usize,
    barcode: &String,
)-> Result<(), Box<dyn Error>>{

    
    // Validate and prepare input data and parameters.
    let (paired_read_file_path_final, read_barcode_file_path_final, single_read_input) = validate_and_assigne_input_reads(r1, r2);
    
    let (mut instrument, mut run) = 
                parse_info_file(
            &find_info_file (&PathBuf::from(info_file), 
                                            &PathBuf::new(), 
                                            &read_barcode_file_path_final
                                            )
                                );
    
    
    if arg_instrument.len() > 0{
        instrument = arg_instrument.clone();
    }
    if arg_run.len() > 0{
        run = arg_run.clone();
    }
    
    let lane: String = match arg_lane.len() > 0{
        true => arg_lane.clone(),
        false => get_lane_from_file(&read_barcode_file_path_final , 1, '_')
    };
    
    
    let (output_directory, report_directory) = prepare_output_report_dir(ouput_dir, report_dir, force);
        
    info!("Output directory: {}", output_directory.display());
    info!("Reports directory: {}", report_directory.display());
    info!("Instrumnet: {}", instrument);
    info!("Run: {}", run);
    info!("Lane: {}", lane);
    info!("Sample Barcode: {}", barcode);
    info!("Compression level: {}. (0 no compression but fast, 12 best compression but slow.)", compression_level);
    

    if illumina_format {
        if lane.len() == 0 || instrument.len() == 0 || run.len() == 0 {
            panic!("Instrument id and run number are required for QC reports and when Illumina format is requested!")
        }
        info!("Read header and Output files: Illumina format.");          
    }else {
        info!("Read header and Output files: MGI format.");
        if lane.len() == 0 {
            panic!("Lane number is required for QC reports!")
        }
    }
    let (sample_name, sb_flowcell, sb_lane, _) = parse_sb_file_name(&read_barcode_file_path_final.file_stem().unwrap().to_string_lossy().to_string());
    let (output_file_r1, output_file_r2) = create_output_file_name(&sample_name, &lane, sample_index, illumina_format);
    
    let barcode_read_output_path = output_directory.join(output_file_r2);
    let paired_read_output_path  = output_directory.join(output_file_r1);
    check_output_file(&barcode_read_output_path, force);
    check_output_file(&paired_read_output_path, force);


    if sb_lane != lane{
        warn!("The lane extracted from the file name ({}) is not the same as the provided lane ({})", sb_lane, lane);
    }
    
    
    info!("Output buffer size: {}", writing_buffer_size);
    if writing_buffer_size < 65536{
        panic!("Writing buffer size '--writing-buffer-size' should be greater than 65536.");
    }
    
    info!("Compression buffer size: {}", compression_buffer_size);
    if compression_buffer_size > writing_buffer_size{
        panic!("Compression buffer size '--compression-buffer-size' should be less than Writing buffer size ('--writing-buffer-size').");
    }
    
    
    let mut compressor = Compressor::new(CompressionLvl::new(compression_level as i32).unwrap());
    
    let reqiured_output_buffer_size = writing_buffer_size + compressor.gzip_compress_bound(compression_buffer_size);

    let relaxed_writing_buffer_size = writing_buffer_size; //( writing_buffer_size as f32 * 0.95 ) as usize;
    // parse sample/index file and get all mismatches
    
    //  Get reads information 

    let mut whole_read_barcode_len: usize = 0;
    let mut whole_paired_read_len: usize = 0;
    let mut tmp_str: String = String::new();
    let mut paired_read_length: usize = 0;
    let only_plus_r1:bool;
    let header_length_r2:usize;
    let mut header_length_r1: usize = 0;

    
    let mut tmp_reader = get_buf_reader(&read_barcode_file_path_final);
    
    tmp_reader.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    header_length_r2 = tmp_str.len();

    let (flowcell, l_position) = get_flowcell_info(&tmp_str);        
    if sb_flowcell != flowcell{
        warn!("The flowcell extracted from the file name ({}) is not the same as the flowcell extracted from the rewad header ({})", sb_flowcell, flowcell);
    }

    let mut sep_position = match memchr(b' ', &tmp_str.as_bytes()) {
        None => { panic!("The header is expected to ends with index information similar to ` *:N:0:********`");},
        Some(loc) => {loc}
    };


    let illumina_header_prefix_str = create_illumina_header_prefix(&instrument, &run, &flowcell);
    let illumina_header_prefix = illumina_header_prefix_str.as_bytes(); 
    
    tmp_str.clear();
    tmp_reader.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    let barcode_read_length: usize = tmp_str.chars().count() - 1;
    
    tmp_str.clear();
    tmp_reader.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    
    let only_plus_r2:bool = tmp_str == "+\n";
    if ! only_plus_r2{
        panic!("Expected read format is not satisified. You can try rerunning using --flexible parameter.");
    }
    tmp_str.clear();
    tmp_reader.read_line(&mut tmp_str).unwrap();
    whole_read_barcode_len += tmp_str.len();
    
    if !single_read_input {
        tmp_str.clear();
        tmp_reader = get_buf_reader(&paired_read_file_path_final);
        tmp_reader.read_line(&mut tmp_str).unwrap();
        whole_paired_read_len += tmp_str.len();
        header_length_r1 = tmp_str.len();

        tmp_str.clear();
        tmp_reader.read_line(&mut tmp_str).unwrap();
        paired_read_length = tmp_str.chars().count() - 1;
        whole_paired_read_len += tmp_str.len();
        
        tmp_str.clear();
        tmp_reader.read_line(&mut tmp_str).unwrap();
        whole_paired_read_len += tmp_str.len();
        only_plus_r1 = tmp_str == "+\n";
        if ! only_plus_r1{
            panic!("Expected read format is not satisified. You can try running demultplex-dynamic command.");
        }
        tmp_str.clear();
        tmp_reader.read_line(&mut tmp_str).unwrap();
        whole_paired_read_len += tmp_str.len();
    }
            
    info!("The length of the read with barcode is: {}", barcode_read_length);
    info!("The length of the paired read is: {}",       paired_read_length);
    
    let mut out_read_barcode_buffer: Vec<u8> = vec![0; reqiured_output_buffer_size];
    let mut out_read_barcode_compression_buffer: Vec<u8> = vec![0; compression_buffer_size];
    let mut out_paired_read_buffer: Vec<u8> = vec![0; reqiured_output_buffer_size];
    let mut out_paired_read_compression_buffer: Vec<u8> = vec![0; compression_buffer_size];
    
   

    let mut read_cntr:u64 = 0;
    
    let shift = if single_read_input{0}else{1};

    let start = Instant::now();
    let dur;
    
    
    

    let mut reader_barcode_read = get_reader(&read_barcode_file_path_final);
    
    let mut reader_paired_read = if !single_read_input {
        let r1 = get_reader(&paired_read_file_path_final);
        Some(r1)
    } else {
        None
    };
    
    let mut buffer_1 = [0; BUFFER_SIZE];  // paired read
    let mut buffer_2 = [0; BUFFER_SIZE]; // read with barcode
    
    let mut read_bytes_1: usize = 0;

    if !single_read_input{
        match reader_paired_read {
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
    let mut actual_sz:usize;
    let mut curr_writer: File;
    
    let mut sample_statistics: Vec<u64> = vec![0; 15];
    let mut out_read_barcode_buffer_last = 0;
    let mut out_paired_read_buffer_last = 0;
    let mut out_read_barcode_compression_buffer_last = 0;
    let mut out_paired_read_compression_buffer_last = 0;
    
    let dynamic_demultiplexing = false;
    let mut curr_umi = vec![b'0'; umi_length];
    let barcode_bytes: &[u8] = barcode.as_bytes();
    let mut curr_mismatch: usize;
    let mut make_warn = true;
    let mut extra_comfort_barcode: usize = whole_read_barcode_len + illumina_header_prefix.len() + 25;
    let mut extra_comfort_paired: usize = whole_paired_read_len + illumina_header_prefix.len() + 25;

    loop {
        //info!("Read: {}", read_cntr);
        read_cntr += 1;
        if dynamic_demultiplexing{
            seq_start = match memchr(b'\n', &buffer_2[header_start..read_bytes_2]) {
                None => { panic!("Something wrong with the input data!");},
                Some(loc) => {header_start + loc + 1}
            };
            sep_position = match memchr(b' ', &buffer_2[header_start..seq_start]) {
                None => { panic!("The header is expected to ends with index information similar to ` *:N:0:********`");},
                Some(loc) => {loc - header_start}
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
            whole_read_barcode_len = 2 * (read_end - header_start + 1);
            extra_comfort_barcode = whole_read_barcode_len + illumina_header_prefix.len() + 25;

        }else{
            seq_start = header_start + header_length_r2;
            plus_start = seq_start + barcode_read_length + 1;
            qual_start = plus_start + 2;
            read_end = barcode_read_length + qual_start;  // \n position.
            if buffer_2[seq_start - 1] != b'\n' || buffer_2[plus_start - 1] != b'\n' ||
                buffer_2[qual_start - 1] != b'\n' || buffer_2[read_end] != b'\n' || buffer_2[sep_position] != b' '{
                panic!("Expected read format is not satisified. You can try rerunning using --flexible parameter.");
            }
            if read_end >= read_bytes_2{
                panic!("Something wrong, read 2 end is beyond the actual buffer size!.");
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
                whole_paired_read_len = 2 * (read_end_pr - header_start_pr + 1);
                extra_comfort_paired = whole_paired_read_len + illumina_header_prefix.len() + 25;

            }else{
                seq_start_pr = header_start_pr + header_length_r1;
                plus_start_pr = seq_start_pr + paired_read_length + 1;
                qual_start_pr = plus_start_pr + 2;
                read_end_pr = paired_read_length + qual_start_pr;
                if buffer_1[seq_start_pr - 1] != b'\n' || buffer_1[plus_start_pr - 1] != b'\n' ||
                    buffer_1[qual_start_pr - 1] != b'\n' || buffer_1[read_end_pr] != b'\n'{
                    panic!("Expected format is not satisified. You can try rerunning using --flexible parameter.");
                }
                if read_end_pr >= read_bytes_1{
                    panic!("Something wrong, read 2 end is beyond the actual buffer size!.");
                }    
            }      
        }
        
        
        /*
            sample_statistics:
            0: r1 count of qs > 30
            1: r2 count os qs > 30
            2: r3 count of (barcode) with qs > 30

            3: r1 count of bases
            4: r2 count of bases
            5: r3 count of bases

            6: qs of r1
            7: qs of r2
            8: qs of r3
            9: read count
            from 10 to whatever mismatches allowed: is the number oif reads with the iteration mismatch
            
        */
        if reporting_level > 0 {
            if barcode_bytes.len() > 0{
                curr_mismatch = hamming_distance(barcode_bytes, &buffer_2[sep_position + header_start + 7..seq_start - 1]).unwrap();
                //println!("{} - {} - {}", curr_mismatch, String::from_utf8(barcode_bytes.to_vec()).unwrap(), 
                //String::from_utf8(buffer_2[sep_position + header_start + 7..seq_start - 1].to_vec()).unwrap());
                if curr_mismatch > 4{
                    curr_mismatch = 4;
                    if make_warn{
                        warn!("This barcode {} (and maybe others) have more than 4 mismatches with the input barocode! Maybe the barcode needs to be reveresed complementary!", 
                        String::from_utf8(buffer_2[sep_position + header_start + 7..seq_start - 1].to_vec()).unwrap());
                        make_warn = false;
                    }  
                }
                sample_statistics[10 + curr_mismatch] += 1; 
                
            }
            
            if ! single_read_input{
                // this is for r1 only if paired end
                for &qs in buffer_1[qual_start_pr..read_end_pr].iter(){
                        if qs >= 63 {
                            sample_statistics[0] += 1;
                            
                        }
                        sample_statistics[6] += qs as u64;
                                            
                }
                
            }

            for &qs in buffer_2[qual_start..read_end - umi_length].iter(){
                if qs >= 63{
                    sample_statistics[shift] += 1;
                }
                sample_statistics[7] += qs as u64;
                
            }

            for &qs in buffer_2[read_end - umi_length..read_end].iter(){
                if qs >= 63 {
                    // this is for r3 or barcode
                    sample_statistics[2] += 1;
                }
                sample_statistics[8] += qs as u64;
                
            }

        }
        
        // writing preperation
        
        // this works for mgi format and unde and ambig and ilumina with a bit of extr
        
        if illumina_format{

            if out_read_barcode_compression_buffer_last + extra_comfort_barcode >= compression_buffer_size {
                actual_sz = compressor.gzip_compress(&out_read_barcode_compression_buffer[..out_read_barcode_compression_buffer_last],
                                    &mut out_read_barcode_buffer[out_read_barcode_buffer_last..]).unwrap();
                out_read_barcode_buffer_last += actual_sz;
    
                out_read_barcode_compression_buffer_last = 0;
    
                    if out_read_barcode_buffer_last >= relaxed_writing_buffer_size {
                        curr_writer = OpenOptions::new()
                                .append(true)
                                .create(true)
                                .open(&barcode_read_output_path)
                                .expect("couldn't create output");
                        curr_writer.write_all(&out_read_barcode_buffer[..out_read_barcode_buffer_last]).unwrap();
                        curr_writer.flush().expect("couldn't flush output");
                        out_read_barcode_buffer_last = 0;
                        
                    }          
            }

            curr_umi[..].copy_from_slice(&buffer_2[plus_start - 1 - umi_length..plus_start - 1]);
            // Illumina format write the heder and skip the header for mgi.
            out_read_barcode_compression_buffer[out_read_barcode_compression_buffer_last..out_read_barcode_compression_buffer_last + illumina_header_prefix.len()].
            copy_from_slice(&illumina_header_prefix[..]);
            out_read_barcode_compression_buffer_last += illumina_header_prefix.len();

            out_read_barcode_compression_buffer_last = write_illumina_header(&mut out_read_barcode_compression_buffer, 
                out_read_barcode_compression_buffer_last, &buffer_2[header_start..seq_start], 
                                    &curr_umi, 
                                    l_position,
                                    sep_position,
                                    true
                                );
        
            header_start = seq_start - 1;

            out_read_barcode_compression_buffer[out_read_barcode_compression_buffer_last..out_read_barcode_compression_buffer_last +  plus_start - header_start - umi_length - 1].
            copy_from_slice(&buffer_2[header_start..plus_start - umi_length - 1]);

            out_read_barcode_compression_buffer_last += plus_start - header_start - umi_length - 1;
                    
            out_read_barcode_compression_buffer[out_read_barcode_compression_buffer_last..out_read_barcode_compression_buffer_last + read_end - plus_start - umi_length + 1].
                        copy_from_slice(&buffer_2[plus_start - 1..read_end - umi_length]);
                
            out_read_barcode_compression_buffer_last += read_end - umi_length - plus_start + 2; 
            out_read_barcode_compression_buffer[out_read_barcode_compression_buffer_last - 1] = b'\n';
                    
            
            if ! single_read_input{
                if out_paired_read_compression_buffer_last + extra_comfort_paired >= compression_buffer_size{                
                    actual_sz = compressor.gzip_compress(&out_paired_read_compression_buffer[.. out_paired_read_compression_buffer_last], 
                        &mut out_paired_read_buffer[out_paired_read_buffer_last..]).unwrap();
                    out_paired_read_buffer_last += actual_sz;
                    out_paired_read_compression_buffer_last = 0;

                    if out_paired_read_buffer_last >= relaxed_writing_buffer_size {
                        curr_writer = OpenOptions::new()
                                .append(true)
                                .create(true)
                                .open(&paired_read_output_path)
                                .expect("couldn't create output");
                        curr_writer.write_all(&out_paired_read_buffer[..out_paired_read_buffer_last]).unwrap();
                        curr_writer.flush().expect("couldn't flush output");
                        out_paired_read_buffer_last = 0;
                    }             
                }

                out_paired_read_compression_buffer[out_paired_read_compression_buffer_last..out_paired_read_compression_buffer_last + illumina_header_prefix.len()].
                copy_from_slice(&illumina_header_prefix[..]);
                out_paired_read_compression_buffer_last += illumina_header_prefix.len();

                out_paired_read_compression_buffer_last = write_illumina_header(&mut out_paired_read_compression_buffer, 
                    out_paired_read_compression_buffer_last, &buffer_1[header_start_pr..seq_start_pr], 
                                        &curr_umi, 
                                        l_position,
                                    sep_position,
                                true);
            
                header_start_pr = seq_start_pr - 1;

                out_paired_read_compression_buffer[out_paired_read_compression_buffer_last..out_paired_read_compression_buffer_last + read_end_pr - header_start_pr + 1].
                copy_from_slice(&buffer_1[header_start_pr..read_end_pr + 1]);
                
                out_paired_read_compression_buffer_last += read_end_pr - header_start_pr + 1;
            
            }        
        
                    
        }

        if ! single_read_input{
            if read_end_pr + whole_paired_read_len >= read_bytes_1 {      
                if read_end_pr < read_bytes_1 - 1 {
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
                        read_bytes(reader, &mut buffer_1, whole_paired_read_len, &mut read_bytes_1);                       
                    },
                    None => panic!("expected sinle end input!")
                };
            }else{
                header_start_pr = read_end_pr + 1;
            }
            
        }

        
        
        if read_end + whole_read_barcode_len >= read_bytes_2{
            if read_end < read_bytes_2 - 1 {
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
            
            read_bytes(&mut reader_barcode_read, &mut buffer_2, whole_read_barcode_len, &mut read_bytes_2);
            
            if read_bytes_2 == 0{
                break;
            }
            
        }else{
            header_start = read_end + 1;
        }
        
        
        
        
    }
    if reporting_level > 0{
        sample_statistics[3] = paired_read_length as u64 * read_cntr;
        sample_statistics[6] -= sample_statistics[3] * 33;
        sample_statistics[shift + 3] = (barcode_read_length as u64 - umi_length as u64) * read_cntr;
        sample_statistics[5] = umi_length as u64 * read_cntr as u64;
        sample_statistics[7] -= sample_statistics[shift + 3] * 33;
        sample_statistics[8] -= sample_statistics[5] * 33;
        sample_statistics[9] = read_cntr;    
    }
    
    if out_read_barcode_compression_buffer_last > 0 {
        actual_sz = compressor.gzip_compress(&out_read_barcode_compression_buffer[..out_read_barcode_compression_buffer_last],
                            &mut out_read_barcode_buffer[out_read_barcode_buffer_last..]).unwrap();
        out_read_barcode_buffer_last += actual_sz;       
    }
    if out_read_barcode_buffer_last > 0 {
        curr_writer = OpenOptions::new()
                .append(true)
                .create(true)
                .open(&barcode_read_output_path)
                .expect("couldn't create output");
        curr_writer.write_all(&out_read_barcode_buffer[..out_read_barcode_buffer_last]).unwrap();
        curr_writer.flush().expect("couldn't flush output");
                
    }   

    if ! single_read_input{
        if out_paired_read_compression_buffer_last > 0{                
            actual_sz = compressor.gzip_compress(&out_paired_read_compression_buffer[.. out_paired_read_compression_buffer_last], 
                &mut out_paired_read_buffer[out_paired_read_buffer_last..]).unwrap();
            out_paired_read_buffer_last += actual_sz;          
        }

        if out_paired_read_buffer_last > 0 {
            curr_writer = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&paired_read_output_path)
                    .expect("couldn't create output");
            curr_writer.write_all(&out_paired_read_buffer[..out_paired_read_buffer_last]).unwrap();
            curr_writer.flush().expect("couldn't flush output");
        }   
    }
    
    dur = start.elapsed();
    
    info!("{} reads were processed in {} secs.", read_cntr, dur.as_secs());
    if reporting_level > 0 {
        for i in (11..15).rev(){
            if sample_statistics[i] == 0{
                sample_statistics.pop();
            }else{
                break;
            }
            
        }
        if barcode.len() == 0{
            sample_statistics[10] = read_cntr;
        }
        let start_logs = Instant::now();
        let mut out_str = String::from(format!("job_number\tsample_id\tr1_qc_30\tr2_qc_30\tr3_qc_30\tr1_bases\tr2_bases\tr3_bases\tr1_qc\tr2_qc\tr3_qc\tall_reads"));
        if sample_statistics.len() > 9 {
            for i in 10..sample_statistics.len(){
                out_str.push('\t');
                out_str.push_str(&(i - 10).to_string());
                out_str.push_str(&"-mismatches".to_string());
            }
        }
        out_str.push('\n');
        out_str.push('.');
        out_str.push('\t');
        out_str.push_str(&sample_name);
        for cnt in 0..sample_statistics.len(){
            out_str.push('\t');
            out_str.push_str(&sample_statistics[cnt].to_string());    
        }
        out_str.push('\n');
        let mut outfile = File::create(report_directory.join(format!("{}.mgikit.sample_stats", sample_name))).expect("couldn't create output");
        outfile.write_all(&out_str.as_bytes()).unwrap();
        let log_dur = start_logs.elapsed();
    
        info!("Writing all logs and reports took {} secs.", log_dur.as_secs());
        
    }
    
    

    Ok(())
}
