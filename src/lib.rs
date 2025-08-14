#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
use clap::ArgMatches;
use crossbeam_channel::{bounded, Receiver, Sender};
use file_utils::{
    get_buf_reader, get_gzip_reader, parallel_reader_decompressor_thread, parallel_reader_thread,
    read_buffers, write_file,
};
use log::{debug, info, warn};
use memchr::{memchr, memchr_iter};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::io::BufRead;
use std::path::Path;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;
use variables::{I5_COLUMN, I7_COLUMN, PROJECT_ID_COLUMN, SAMPLE_COLUMN};

// my modules
mod file_utils;
mod formater;
mod hardware_resources;
mod report_manager;
mod run_manager;
mod sample_data;
mod sample_manager;
mod variables;

pub use crate::hardware_resources::{get_available_memory, get_cpus};
pub use crate::sample_data::*;
pub use formater::{parse_sb_file_name, ReformatedSample};
pub use report_manager::ReportManager;
pub use run_manager::RunManager;
pub use sample_manager::{get_all_mismatches, reverse_complement, SampleManager};

const BUFFER_SIZE: usize = 1 << 22;
const RAW_BUFFER_SIZE: usize = 1 << 23;

fn hamming_distance(bytes1: &[u8], bytes2: &[u8]) -> Result<usize, &'static str> {
    if bytes1.len() != bytes2.len() {
        return Err("Byte slices must be of equal length");
    }

    Ok(bytes1
        .iter()
        .zip(bytes2.iter())
        .filter(|(b1, b2)| b1 != b2)
        .count())
}

pub fn add_match(
    curr_template: String,
    sample_itr: usize,
    matches_stat: &mut HashMap<String, Vec<usize>>,
    sample_list_ln: usize,
) {
    match matches_stat.get_mut(&curr_template) {
        Some(sample_reads) => {
            sample_reads[sample_itr] += 1;
        }
        None => {
            let mut tmp = vec![0; sample_list_ln];
            tmp[sample_itr] += 1;
            matches_stat.insert(curr_template, tmp);
        }
    };
}

fn find_matches(
    sample_itr: usize,
    matches_stat: &mut HashMap<String, Vec<usize>>,
    sample_list_ln: usize,
    read_barcode_seq: &String,
    sample_first_indx: &String,
    sample_first_indx_rc: &String,
    sample_second_indx: &String,
    sample_second_indx_rc: &String,
) {
    //println!("{} - {} - {} - {}  -> {}", sample_first_indx, sample_first_indx_rc, sample_second_indx, sample_second_indx_rc, read_barcode_seq);
    let indexes = match sample_second_indx.len() < 2 {
        true => vec![sample_first_indx, sample_first_indx_rc],
        false => vec![
            sample_first_indx,
            sample_first_indx_rc,
            sample_second_indx,
            sample_second_indx_rc,
        ],
    };
    //println!("{:?}", indexes);
    for i in 0..2 {
        let first_matches: Vec<usize> = read_barcode_seq
            .match_indices(indexes[i])
            .map(|(i, _)| i)
            .collect();
        //println!("---------------------");
        //println!("{} -> {}, matches first: {:?}", i, indexes[i], first_matches);

        if first_matches.len() > 0 {
            for j in 2..indexes.len() {
                let second_matches: Vec<usize> = read_barcode_seq
                    .match_indices(indexes[j])
                    .map(|(i, _)| i)
                    .collect();
                //println!("{} -> {}, matches second: {:?}", j, indexes[j], second_matches);
                for second_match in &second_matches {
                    for first_match in &first_matches {
                        //println!("sample - {}: comparing  {}  ->  {}    {}:{}",sample_itr, first_match, second_match, i, j);

                        if *second_match >= first_match + sample_first_indx.len()
                            || *first_match >= second_match + sample_second_indx.len()
                        {
                            add_match(
                                format!("{}:{}_{}:{}", first_match, second_match, i, j - 2),
                                sample_itr,
                                matches_stat,
                                sample_list_ln,
                            );
                        }
                    }
                }
            }

            if indexes.len() == 2 {
                for first_match in first_matches {
                    add_match(
                        format!("{}_{}:.", first_match, i),
                        sample_itr,
                        matches_stat,
                        sample_list_ln,
                    );
                }
            }
        }
    }
}

pub fn get_mgikit_template(
    initial_template: &String,
    umi: bool,
    barcode_length: usize,
    i7_len: usize,
    i5_len: usize,
) -> (String, String, String) {
    let mut final_template: Vec<String> = Vec::new();
    let template_elements: Vec<String> = initial_template.split('_').map(String::from).collect();
    let rc_ls: Vec<String> = template_elements[1].split(':').map(String::from).collect();
    let mut index_ls: Vec<usize> = template_elements[0]
        .split(':')
        .map(|it| it.parse().unwrap())
        .collect();

    let mut shift = 0;
    let no_swap: bool = index_ls.len() == 1 || index_ls[0] < index_ls[1];
    index_ls.sort();

    let mut umi_loc = 0;
    let mut umi_size = 0;
    let mut item_added = 0;
    //println!("{:?}", index_ls);
    if index_ls[0] > 0 {
        final_template.push(format!("--{}", index_ls[0]));
        shift += index_ls[0];
        umi_size = index_ls[0];
        item_added += 1;
    }
    //println!("1-   {:?}", final_template);
    if no_swap {
        final_template.push(format!("i7{}", i7_len));
        shift += i7_len;
    } else {
        final_template.push(format!("i5{}", i5_len));
        shift += i5_len;
    }
    item_added += 1;
    //println!("2-   {:?}", final_template);

    if index_ls.len() == 2 {
        if shift < index_ls[1] {
            final_template.push(format!("--{}", index_ls[1] - shift));
            if index_ls[1] - shift > umi_size {
                umi_loc = item_added;
                umi_size = index_ls[1] - shift;
            }
            shift = index_ls[1];
            item_added += 1;
        }

        if no_swap {
            final_template.push(format!("i5{}", i5_len));
            shift += i5_len;
        } else {
            final_template.push(format!("i7{}", i7_len));
            shift += i7_len;
        }
        item_added += 1;
    }
    //println!("3-   {:?}", final_template);

    if shift < barcode_length {
        final_template.push(format!("--{}", barcode_length - shift));
        if barcode_length - shift > umi_size {
            umi_loc = item_added;
        }
    }
    //println!("4-   {:?}", final_template);

    if umi {
        //println!("{}  ->  {}", final_template[umi_loc], umi_loc);
        final_template[umi_loc] = final_template[umi_loc].replace("--", "um");
    }
    //println!("5-   {:?}", final_template);

    (
        final_template.join(":"),
        rc_ls[0].to_owned(),
        rc_ls[1].to_owned(),
    )
}

pub fn find_matching_sample(
    all_template_data: &Vec<(
        u32,
        HashSet<String>,
        HashSet<String>,
        String,
        HashMap<String, (usize, HashMap<String, usize>)>,
        bool,
        [usize; 10],
    )>,
    mismatches_dic_i7: &Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>,
    mismatches_dic_i5: &Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>,
    read_barcode: &[u8],
    allowed_mismatches: usize,
    all_index_allowed_mismatches: usize,
    total_samples: usize,
    comprehensive_scan: bool,
    i7_rc: bool,
) -> (usize, usize, String, String) {
    let mut template_itr = 0;
    let undetermined_label_id = total_samples - 2;
    let ambiguous_label_id = total_samples - 1;
    let mut sample_id: usize = undetermined_label_id;
    let mut curr_umi: String = String::new();
    let mut curr_barcode = String::new();
    let mut latest_mismatch = usize::MAX;
    let mut curr_mismatch = usize::MAX;

    for template_details in all_template_data {
        let sample_info = &template_details.4;
        let indexes_info = template_details.6;

        let all_i7_mismatches = &mismatches_dic_i7[template_itr];
        let all_i5_mismatches = &mismatches_dic_i5[template_itr];

        //println!("{:?}", indexes_info);
        //println!("{}", template_details.5);
        let extract_umi = if indexes_info[6] == 1 { true } else { false };

        match all_i7_mismatches.get(
            &read_barcode[read_barcode.len() - indexes_info[2]
                ..read_barcode.len() - indexes_info[2] + indexes_info[1]],
        ) {
            Some(i7_matches) => {
                //info!("{:?}", i7_matches);
                if i7_matches.1 <= allowed_mismatches {
                    if template_details.5 {
                        match all_i5_mismatches.get(
                            &read_barcode[read_barcode.len() - indexes_info[5]
                                ..read_barcode.len() - indexes_info[5] + indexes_info[4]],
                        ) {
                            Some(i5_matches) => {
                                //info!("{:?} and {:?}", i7_matches, i5_matches);
                                for i7_match_itr in 0..i7_matches.0.len() {
                                    for i5_match_itr in 0..i5_matches.0.len() {
                                        if i5_matches.1 > allowed_mismatches
                                            || latest_mismatch < i7_matches.1 + i5_matches.1
                                        {
                                            continue;
                                        }
                                        curr_mismatch = i7_matches.1 + i5_matches.1;
                                        if curr_mismatch <= all_index_allowed_mismatches {
                                            //println!("{}", i7_matches.0[i7_match_itr]);
                                            match sample_info
                                                .get(i7_matches.0[i7_match_itr])
                                                .unwrap()
                                                .1
                                                .get(i5_matches.0[i5_match_itr])
                                            {
                                                Some(i5_info) => {
                                                    if sample_id >= total_samples
                                                        || latest_mismatch > curr_mismatch
                                                    {
                                                        sample_id = *i5_info;
                                                        latest_mismatch = curr_mismatch;
                                                        curr_barcode = unsafe {
                                                            String::from_utf8_unchecked(
                                                                read_barcode[read_barcode.len()
                                                                    - indexes_info[2]
                                                                    ..read_barcode.len()
                                                                        - indexes_info[2]
                                                                        + indexes_info[1]]
                                                                    .to_vec(),
                                                            )
                                                        };
                                                        curr_barcode.push('+');
                                                        curr_barcode.push_str(
                                                            &(unsafe {
                                                                String::from_utf8_unchecked(
                                                                    read_barcode[read_barcode.len()
                                                                        - indexes_info[5]
                                                                        ..read_barcode.len()
                                                                            - indexes_info[5]
                                                                            + indexes_info[4]]
                                                                        .to_vec(),
                                                                )
                                                            }),
                                                        );

                                                        if extract_umi {
                                                            curr_umi = String::from(":");
                                                            if i7_rc {
                                                                curr_umi.push_str(
                                                                    &reverse_complement(
                                                                        &(unsafe {
                                                                            String::from_utf8_unchecked(
                                                                                read_barcode[
                                                                                    read_barcode.len() -
                                                                                        indexes_info
                                                                                            [8]..read_barcode.len() -
                                                                                        indexes_info
                                                                                            [8] +
                                                                                        indexes_info
                                                                                            [7]
                                                                                ].to_vec()
                                                                            )
                                                                        })
                                                                    ).unwrap()
                                                                );
                                                            } else {
                                                                curr_umi.push_str(
                                                                    &(unsafe {
                                                                        String::from_utf8_unchecked(
                                                                            read_barcode[
                                                                                read_barcode.len() -
                                                                                    indexes_info
                                                                                        [8]..read_barcode.len() -
                                                                                    indexes_info
                                                                                        [8] +
                                                                                    indexes_info[7]
                                                                            ].to_vec()
                                                                        )
                                                                    })
                                                                );
                                                            }
                                                        }
                                                    } else {
                                                        sample_id = ambiguous_label_id;
                                                        break;
                                                    }
                                                }
                                                None => {
                                                    //break
                                                }
                                            };
                                        }
                                    }
                                }
                            }
                            None => {}
                        }
                    } else {
                        if latest_mismatch < i7_matches.1 {
                            continue;
                        } else if latest_mismatch > i7_matches.1 && i7_matches.0.len() < 2 {
                            curr_mismatch = i7_matches.1;
                            latest_mismatch = curr_mismatch;
                            sample_id = match template_details.4.get(i7_matches.0[0]) {
                                Some(tmp) => tmp.0,
                                None => {
                                    panic!("There should be a value here");
                                }
                            };
                            //barcode_length = indexes_info[9];
                            curr_barcode = unsafe {
                                String::from_utf8_unchecked(
                                    read_barcode[read_barcode.len() - indexes_info[2]
                                        ..read_barcode.len() - indexes_info[2] + indexes_info[1]]
                                        .to_vec(),
                                )
                            };
                            if extract_umi {
                                curr_umi = String::from(":");
                                curr_umi.push_str(
                                    &(unsafe {
                                        String::from_utf8_unchecked(
                                            read_barcode[read_barcode.len() - indexes_info[8]
                                                ..read_barcode.len() - indexes_info[8]
                                                    + indexes_info[7]]
                                                .to_vec(),
                                        )
                                    }),
                                );
                                /*if i7_rc {
                                    curr_umi.push_str(
                                        &reverse_complement(
                                            &(unsafe {
                                                String::from_utf8_unchecked(
                                                    read_barcode[
                                                        read_barcode.len() -
                                                            indexes_info
                                                                [8]..read_barcode.len() -
                                                            indexes_info[8] +
                                                            indexes_info[7]
                                                    ].to_vec()
                                                )
                                            })
                                        ).unwrap()
                                    );
                                }*/
                            }
                        } else {
                            sample_id = ambiguous_label_id;
                        }
                    }
                }
            }
            None => {}
        }

        if sample_id != undetermined_label_id && sample_id < total_samples && !comprehensive_scan {
            break;
        }

        template_itr += 1;
    }
    if sample_id >= undetermined_label_id {
        curr_mismatch = 0;
    }
    (sample_id, curr_mismatch, curr_barcode, curr_umi)
}

fn copy_within_a_slice<T: Clone>(v: &mut [T], from: usize, to: usize, len: usize) {
    if len == 0 {
        return;
    }
    if from > to {
        let (dst, src) = v.split_at_mut(from);
        dst[to..to + len].clone_from_slice(&src[..len]);
    } else {
        let (src, dst) = v.split_at_mut(to);
        dst[..len].clone_from_slice(&src[from..from + len]);
    }
}

fn sum_qc(qc_seq: &[u8], check_content: bool) -> (u64, u64) {
    let mut qc_total: u64 = 0;
    let mut high_qc: u64 = 0;
    for &qs in qc_seq.iter() {
        if check_content {
            if qs < 33 || qs > 73 {
                panic!("Reverse read quality scores should be between 33 and 73. the detected quality score is {}", qs);
            }
        }
        if qs >= 63 {
            high_qc += 1;
        }
        qc_total += qs as u64;
    }
    //debug!("quality res: {} - {} - {}", qc_seq.len(), qc_total, high_qc);
    (qc_total, high_qc)
}

pub fn check_read_content(read_seq: &[u8], seq_start: usize, plus_start: usize, qc_start: usize) {
    if read_seq[0] != b'@' {
        panic!(
            "Read header must starts with '@'! The detected header is : {}",
            String::from_utf8(read_seq.to_vec()).unwrap()
        );
    }
    for &nec in read_seq[seq_start..plus_start - 1].iter() {
        if nec != b'A'
            && nec != b'C'
            && nec != b'G'
            && nec != b'T'
            && nec != b'N'
            && nec != b'a'
            && nec != b'c'
            && nec != b'g'
            && nec != b't'
            && nec != b'n'
        {
            panic!(
                "Seqeunce bases need to be in [A, C, G, T, a, c, g, t, N, n]! Found {}.",
                nec
            );
        }
    }
    for &qs in read_seq[qc_start..read_seq.len()].iter() {
        if qs > 73 || qs < 33 {
            panic!(
                "Quality score must be between [0 and 40]! Found {}!",
                qs - 33
            );
        }
    }
}

fn demultiplex(
    sample_manager: &SampleManager,
    run_manager: &RunManager,
    buffer_info: &BufferInfo,
    allowed_mismatches: usize,
    reporting_level: usize,
    all_index_error: bool,
    reader_threads: usize,
    processing_threads: usize,
) -> Result<ReportManager, Box<dyn Error>> {
    let start = Instant::now();
    //let dur;

    info!(
        "Comprehensive scan mode: {}",
        run_manager.comprehensive_scan()
    );
    info!(
        "Allowed mismatches when finding the index are: {}",
        allowed_mismatches
    );

    if all_index_error {
        info!(
            "The allowed mismatches will be compared to the total mismatches for all indicies combined."
        );
    } else {
        info!("The allowed mismatches will be considered per index.");
    }
    if allowed_mismatches > 2 {
        warn!(
            "It is not recommended to allow more than 2 mismatches per index! This might decrease the accuracy of the demultipexing!"
        );
    }

    info!("Reads that match with multiple samples will be saved in ambiguous_read1/2.fastq file");
    info!("The paired read files are assumed to contain the paired reads in the same order!");

    let all_template_data = sample_manager.all_template_data();
    let total_samples = sample_manager.get_sample_count();

    if run_manager.force() {
        clean_output_directory(
            sample_manager,
            run_manager,
            buffer_info,
            run_manager.read2_has_sequence(),
            run_manager.illumina_format(),
        );
    }

    let (full_sender_rb, full_receiver_rb) = bounded(processing_threads * 2);
    let (empty_sender_rb, empty_receiver_rb) = bounded(processing_threads * 2);
    let (full_sender_rp, full_receiver_rp) = bounded(processing_threads * 2);
    let (empty_sender_rp, empty_receiver_rp) = bounded(processing_threads * 2);
    let (full_sender, full_receiver) = bounded(processing_threads * 2);

    /*
    empty_receiver_rb/empty_receiver_rp
    full_sender_rb/full_sender_rp
    full_receiver_rb/full_receiver_rp
    full_sender


    workers
    full_receiver
    empty_sender_rb
    empty_sender_rp


    empty_raw_receiver_rp/empty_raw_receiver_rb
    full_raw_sender_rp/full_raw_sender_rb

    decoder
    full_raw_receiver_rp/full_raw_receiver_rb
    empty_receiver_rb/empty_receiver_rp
    full_sender_rb/full_sender_rp
    empty_raw_sender_rb/empty_raw_sender_rp
    full_receiver_rb/full_receiver_rp
    full_sender
    */

    let (full_raw_sender_rb, full_raw_receiver_rb) = bounded(processing_threads * 2);
    let (empty_raw_sender_rb, empty_raw_receiver_rb) = bounded(processing_threads * 2);
    let (full_raw_sender_rp, full_raw_receiver_rp) = bounded(processing_threads * 2);
    let (empty_raw_sender_rp, empty_raw_receiver_rp) = bounded(processing_threads * 2);

    let mut barcode_process_master = !run_manager.paired_read_input() || reader_threads == 1;

    let batch_size: usize = (if run_manager.barcode_read_len() > run_manager.paired_read_len() {
        ((BUFFER_SIZE as f32) / (run_manager.barcode_read_len() as f32)) * 0.95
    } else {
        barcode_process_master = true;
        ((BUFFER_SIZE as f32) / (run_manager.paired_read_len() as f32)) * 0.95
    }) as usize;
    debug!(
        "Barceode data reader is master reader: {}",
        barcode_process_master
    );
    let reader_handler = if reader_threads > 0 {
        info!(
            "Parallel readers, processing batches of {} reads.",
            batch_size
        );
        for _ in 0..processing_threads * 2 {
            match empty_sender_rb.send((0, vec![0u8; BUFFER_SIZE])) {
                Ok(_) => {}
                Err(e) => println!("0- Error Sending: {}", e),
            }
            match empty_raw_sender_rb.send((0, vec![0u8; RAW_BUFFER_SIZE])) {
                Ok(_) => {}
                Err(e) => println!("0- Error Sending: {}", e),
            }
            if run_manager.paired_read_input() {
                match empty_sender_rp.send((0, vec![0u8; BUFFER_SIZE])) {
                    Ok(_) => {}
                    Err(e) => println!("0- Error Sending: {}", e),
                }
                match empty_raw_sender_rp.send((0, vec![0u8; RAW_BUFFER_SIZE])) {
                    Ok(_) => {}
                    Err(e) => println!("0- Error Sending: {}", e),
                }
            }
        }
        if (reader_threads == 4 && run_manager.paired_read_input())
            || (reader_threads == 2 && !run_manager.paired_read_input())
        {
            Some(parallel_reader_decompressor_thread(
                run_manager.barcode_reads().clone(),
                batch_size,
                empty_raw_sender_rb,
                empty_raw_receiver_rb,
                full_raw_sender_rb,
                full_raw_receiver_rb,
                full_receiver_rb.clone(),
                full_receiver_rp.clone(),
                full_sender_rb.clone(),
                empty_receiver_rb.clone(),
                empty_sender_rb.clone(),
                full_sender.clone(),
                processing_threads,
                BUFFER_SIZE,
                run_manager.paired_read_input(),
                barcode_process_master,
            ))
        } else {
            Some(parallel_reader_thread(
                run_manager.paired_reads().clone(),
                run_manager.barcode_reads().clone(),
                batch_size,
                true,
                reader_threads == 1 && run_manager.paired_read_input(),
                full_sender_rb.clone(),
                full_sender_rp.clone(),
                empty_sender_rb.clone(),
                empty_sender_rp.clone(),
                full_receiver_rb.clone(),
                full_receiver_rp.clone(),
                empty_receiver_rb.clone(),
                empty_receiver_rp.clone(),
                full_sender.clone(),
                processing_threads,
                BUFFER_SIZE,
                run_manager.paired_read_input(),
                barcode_process_master,
            ))
        }
    } else {
        None
    };

    let reader_handler_secondary = if run_manager.paired_read_input() && reader_threads > 1 {
        if reader_threads == 4 {
            Some(parallel_reader_decompressor_thread(
                run_manager.paired_reads().clone(),
                batch_size,
                empty_raw_sender_rp,
                empty_raw_receiver_rp,
                full_raw_sender_rp,
                full_raw_receiver_rp,
                full_receiver_rb.clone(),
                full_receiver_rp.clone(),
                full_sender_rp.clone(),
                empty_receiver_rp.clone(),
                empty_sender_rp.clone(),
                full_sender.clone(),
                processing_threads,
                BUFFER_SIZE,
                true,
                !barcode_process_master,
            ))
        } else {
            Some(parallel_reader_thread(
                run_manager.paired_reads().clone(),
                run_manager.barcode_reads().clone(),
                batch_size,
                false,
                true,
                full_sender_rb.clone(),
                full_sender_rp.clone(),
                empty_sender_rb.clone(),
                empty_sender_rp.clone(),
                full_receiver_rb.clone(),
                full_receiver_rp.clone(),
                empty_receiver_rb.clone(),
                empty_receiver_rp.clone(),
                full_sender.clone(),
                processing_threads,
                BUFFER_SIZE,
                true,
                !barcode_process_master,
            ))
        }
    } else {
        None
    };

    let report_manager_arc = Arc::new(Mutex::new(ReportManager::new(
        total_samples,
        allowed_mismatches,
    )));
    let mut processor_pool = Vec::new();
    let mut sample_lockes = Vec::new();
    for _ in 0..total_samples {
        sample_lockes.push(Mutex::new(false));
    }
    let samples_locks_arc = Arc::new(sample_lockes);

    for th_id in 0..processing_threads - 1 {
        let thread_name = format!("worker-{}", th_id);
        let full_receiver = full_receiver.clone();
        let empty_sender_rb = empty_sender_rb.clone();
        let empty_sender_rp = empty_sender_rp.clone();
        let sample_manager = sample_manager.clone();
        let run_manager = run_manager.clone();
        let buffer_info = buffer_info.clone();
        let all_template_data = all_template_data.clone();
        let report_manager_arc = report_manager_arc.clone();
        let samples_locks_arc: Arc<Vec<Mutex<bool>>> = samples_locks_arc.clone();
        processor_pool.push(
            thread::Builder::new()
                .name(thread_name.clone())
                .spawn(move || {
                    let curr_report_manager = analyse_fastq(
                        &sample_manager,
                        &run_manager,
                        &buffer_info,
                        all_template_data.clone(),
                        reporting_level,
                        allowed_mismatches,
                        all_index_error,
                        BUFFER_SIZE,
                        full_receiver,
                        empty_sender_rb,
                        empty_sender_rp,
                        reader_threads > 0,
                        samples_locks_arc,
                        true,
                        &ReformatedSample::default(),
                        vec![false; 10],
                    );

                    let mut report_manager = report_manager_arc.lock().unwrap();
                    report_manager.update(&curr_report_manager);
                }),
        );
    }

    let curr_report_manager = analyse_fastq(
        &sample_manager,
        &run_manager,
        &buffer_info,
        all_template_data.clone(),
        reporting_level,
        allowed_mismatches,
        all_index_error,
        BUFFER_SIZE,
        full_receiver,
        empty_sender_rb,
        empty_sender_rp,
        reader_threads > 0,
        samples_locks_arc,
        true,
        &ReformatedSample::default(),
        vec![true; 10],
    );

    if reader_threads > 0 {
        let _ = match reader_handler {
            Some(handler) => handler.join().unwrap(),
            None => {
                panic!("there must be a reader thread here!")
            }
        };
        if reader_threads == 2 {
            let _ = match reader_handler_secondary {
                Some(handler) => handler.join().unwrap(),
                None => {
                    panic!("there must be a secondary reader thread here!")
                }
            };
        }

        for handler in processor_pool {
            handler.unwrap().join().unwrap();
        }
    }

    let mut report_manager = match Arc::try_unwrap(report_manager_arc) {
        Ok(mutex) => match mutex.into_inner() {
            Ok(manager) => manager,
            Err(_) => panic!("Report Manager Mutex poisoned"),
        },
        Err(_) => panic!("Report Manager Arc still has multiple owners"),
    };

    report_manager.update(&curr_report_manager);

    let dur = start.elapsed();

    info!(
        "{} reads were processed in {} secs.",
        report_manager.get_total_reads(),
        dur.as_secs()
    );

    Ok(report_manager)
}

fn process_buffer(
    run_manager: &RunManager,
    _buffer_info: &BufferInfo,
    report_manager: &mut ReportManager,
    all_index_allowed_mismatches: usize,
    sample_manager: &SampleManager,
    all_template_data: &Vec<(
        u32,
        HashSet<String>,
        HashSet<String>,
        String,
        HashMap<String, (usize, HashMap<String, usize>)>,
        bool,
        [usize; 10],
    )>,
    samples_reads: &mut Vec<SampleData>,
    reporting_level: usize,
    allowed_mismatches: usize,
    buffer_1: &[u8],
    buffer_2: &[u8],
    mismatches_dic_i7: &Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>,
    mismatches_dic_i5: &Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>>,
    samples_locks: &Arc<Vec<Mutex<bool>>>,
    lines_rb: Vec<usize>,
    lines_rp: Vec<usize>,
    demultiplex: bool,
    reformated_sample: &ReformatedSample,
    warning_ls: &mut Vec<bool>,
) -> (usize, usize) {
    let l_position: usize = run_manager.l_position();
    let total_samples: usize = sample_manager.get_sample_count();
    let undetermined_label_id = total_samples - 2;
    let ambiguous_label_id = total_samples - 1;
    let shift: usize = if run_manager.paired_read_input() {
        1
    } else {
        0
    };
    let barcode_length: usize = if demultiplex {
        all_template_data[0].6[9]
    } else {
        reformated_sample.umi_length()
    };
    let writen_barcode_length: usize = match run_manager.keep_barcode() {
        false => barcode_length,
        true => 0,
    };
    let writing_samples: &Vec<usize> = sample_manager.writing_samples();
    //let mut curr_template;
    let mut header_start = 0;
    let mut sample_id;
    let mut barcode_read_illumina_header_start: usize;
    let mut curr_mismatch: usize;
    let mut curr_umi = String::new();
    let mut curr_barcode;
    let mut read_end: usize;
    let mut read_end_pr: usize = 0;
    let mut seq_start: usize;
    let mut plus_start: usize;
    let mut qual_start: usize;
    let mut header_start_pr: usize = 0;
    let mut seq_start_pr: usize = 0;
    let mut plus_start_pr: usize = 0;
    let mut qual_start_pr: usize = 0;

    let mut curr_writing_sample: usize;
    let undertmined_threshold_check = 5000;
    let mut reached_an_end = false;
    let mut i7_rc;
    let mut read_cntr: u64 = 0;
    let check_content = run_manager.check_content();

    //let mut lines_rb: Memchr = memchr_iter(b'\n', &buffer_2);
    //let mut lines_rp: Memchr = memchr_iter(b'\n', &buffer_1);
    //debug!("lines rb: {}    lines rp:  {} ", lines_rb.count(), lines_rp.count());

    let illumina_format = run_manager.illumina_format();
    let comprehensive_scan = run_manager.comprehensive_scan();
    let mut lines_rp = lines_rp.into_iter();
    let mut lines_rb = lines_rb.into_iter();
    let mut sep_position: usize;
    let mut header_shift;
    let mut tail_offset;
    let mut header_info: String;
    let mut raw_shift: usize;
    loop {
        header_shift = 0;
        raw_shift = 0;
        //info!("Read: {}", read_cntr);
        seq_start = match lines_rb.next() {
            Some(loc) => loc + 1,
            None => {
                reached_an_end = true;
                0
            }
        };
        plus_start = match lines_rb.next() {
            Some(loc) => loc + 1,
            None => {
                reached_an_end = true;
                0
            }
        };
        qual_start = match lines_rb.next() {
            Some(loc) => loc + 1,
            None => {
                reached_an_end = true;
                0
            }
        };
        read_end = match lines_rb.next() {
            Some(loc) => loc,
            None => {
                reached_an_end = true;
                0
            }
        };

        //debug!("barcode read: {}, start: {}, {}, {}, {}.", read_cntr, header_start, seq_start, plus_start, read_end);
        if run_manager.paired_read_input() {
            seq_start_pr = match lines_rp.next() {
                Some(loc) => loc + 1,
                None => {
                    reached_an_end = true;
                    0
                }
            };
            plus_start_pr = match lines_rp.next() {
                Some(loc) => loc + 1,
                None => {
                    reached_an_end = true;
                    0
                }
            };
            qual_start_pr = match lines_rp.next() {
                Some(loc) => loc + 1,
                None => {
                    reached_an_end = true;
                    0
                }
            };
            read_end_pr = match lines_rp.next() {
                Some(loc) => loc,
                None => {
                    reached_an_end = true;
                    0
                }
            };
            //debug!("paired read: {}, start: {}, {}, {}, {}.", read_cntr, header_start_pr, seq_start_pr, plus_start_pr, read_end_pr);
        }

        if reached_an_end {
            debug!(
                "reached end - buffer 2 - bytes: {}  - header start: {}",
                buffer_2.len(),
                header_start
            );
            debug!(
                "reached end - buffer 1 - bytes: {}  - header start: {}",
                buffer_1.len(),
                header_start_pr
            );
            break;
        } else {
            i7_rc = false;

            if demultiplex {
                sep_position = seq_start - header_start - 3;
                (sample_id, curr_mismatch, curr_barcode, curr_umi) = find_matching_sample(
                    &all_template_data,
                    &mismatches_dic_i7,
                    &mismatches_dic_i5,
                    &buffer_2[plus_start - 1 - barcode_length..plus_start - 1],
                    allowed_mismatches,
                    all_index_allowed_mismatches,
                    total_samples,
                    comprehensive_scan,
                    i7_rc,
                );
                tail_offset = curr_barcode.len() + 6;
                if sample_id >= total_samples {
                    sample_id = undetermined_label_id;
                }

                if sample_id == undetermined_label_id || sample_id == ambiguous_label_id {
                    if warning_ls[0]
                        && read_cntr > undertmined_threshold_check
                        && report_manager.get_sample_reads(sample_id)
                            > ((0.75 * (read_cntr as f64)) as u64)
                    {
                        if run_manager.ignore_undetermined() {
                            warn!(
                                "{}\nAll reads: {}, Undetermined reads: {}",
                                "Seems that there is an issue with the input. Most of the reads are undetermined!",
                                read_cntr,
                                report_manager.get_sample_reads(sample_id)
                            );
                            warning_ls[0] = false;
                        } else {
                            panic!(
                                "{}\nAll reads: {}, Undetermined reads: {}",
                                "Seems that there is an issue with the input. Most of the reads are undetermined!",
                                read_cntr,
                                report_manager.get_sample_reads(sample_id)
                            );
                        }
                    }

                    if all_template_data.len() == 1 {
                        let indexes_info = all_template_data[0].6;
                        curr_barcode = unsafe {
                            String::from_utf8_unchecked(
                                buffer_2[plus_start - indexes_info[2] - 1
                                    ..plus_start - indexes_info[2] + indexes_info[1] - 1]
                                    .to_vec(),
                            )
                        };

                        if indexes_info[3] == 1 {
                            // i5 exist
                            curr_barcode.push('+');
                            curr_barcode.push_str(
                                &(unsafe {
                                    String::from_utf8_unchecked(
                                        buffer_2[plus_start - indexes_info[5] - 1
                                            ..plus_start - indexes_info[5] + indexes_info[4] - 1]
                                            .to_vec(),
                                    )
                                }),
                            );
                        }
                    } else {
                        curr_barcode = unsafe {
                            String::from_utf8_unchecked(
                                buffer_2[plus_start - barcode_length - 1..plus_start - 1].to_vec(),
                            )
                        };
                    }
                } else if curr_barcode.len() == 0 && run_manager.mgi_data() {
                    panic!("Barcode must be captured here!");
                }

                if reporting_level > 1 {
                    if sample_id == undetermined_label_id {
                        report_manager.update_undetermined(curr_barcode.clone(), 1);
                    } else if sample_id == ambiguous_label_id {
                        report_manager.update_ambiguous(curr_barcode.clone(), 1);
                    }
                }
            } else {
                sep_position = match memchr(b' ', &buffer_2[header_start..seq_start]) {
                    None => {
                        if warning_ls[1] {
                            warn!("You might be using MGI raw data");
                            warning_ls[1] = false;
                        }

                        //panic!("The header is expected to end with index information similar to ` *:N:0:********`");
                        raw_shift = 2;
                        seq_start - header_start - 3
                    }
                    Some(loc) => loc,
                };
                curr_mismatch = if reformated_sample.barcode().len() == 0
                    || sep_position == seq_start - header_start - 3
                {
                    0
                } else {
                    hamming_distance(
                        reformated_sample.barcode().as_bytes(),
                        &buffer_2[sep_position + header_start + 7..seq_start - 1],
                    )
                    .unwrap()
                };
                sample_id = 0; //reformated_sample.sample_index();

                //curr_barcode = "";
                if reformated_sample.umi_length() > 0 {
                    unsafe {
                        curr_umi = String::from_utf8_unchecked(
                            buffer_2
                                [plus_start - 1 - reformated_sample.umi_length()..plus_start - 1]
                                .to_vec(),
                        );
                    }
                }

                curr_barcode = curr_umi.clone();
                //unsafe {String::from_utf8_unchecked(buffer_2[sep_position..seq_start - 1].to_vec()).split(":").last().unwrap().to_string();};
                header_shift = seq_start - sep_position - header_start - 3;
                tail_offset = header_shift + 1;

                if reformated_sample.barcode().len() > 0
                    && sep_position == seq_start - header_start - 3
                {
                    header_shift = 0;
                    tail_offset = reformated_sample.barcode().len() + 6;
                }
                //println!("{} - {} - {}", curr_mismatch, String::from_utf8(barcode_bytes.to_vec()).unwrap(),
                //String::from_utf8(buffer_2[sep_position + header_start + 7..seq_start - 1].to_vec()).unwrap());
                if curr_mismatch > 4 {
                    curr_mismatch = 4;
                    if warning_ls[2] {
                        warn!(
                            "This barcode {} (and maybe others) have more than 4 mismatches with the input barocode! Maybe the barcode needs to be reveresed complementary!",
                            String::from_utf8(
                                buffer_2[sep_position + header_start + 7..seq_start - 1].to_vec()
                            ).unwrap()
                        );
                        warning_ls[2] = false;
                    }
                }
            }

            if check_content {
                if run_manager.paired_read_input() {
                    check_read_content(
                        &buffer_1[header_start_pr..read_end_pr],
                        seq_start_pr - header_start_pr,
                        plus_start_pr - header_start_pr,
                        qual_start_pr - header_start_pr,
                    );

                    if run_manager.mgi_data()
                        && hamming_distance(
                            &buffer_1[header_start_pr..seq_start_pr - 2],
                            &buffer_2[header_start..seq_start - 2],
                        )
                        .unwrap()
                            > 0
                    {
                        panic!(
                            "Headers seem to be in differnet orders: {}  -  {}",
                            String::from_utf8(buffer_1[header_start_pr..seq_start_pr - 2].to_vec())
                                .unwrap(),
                            String::from_utf8(buffer_2[header_start..seq_start - 2].to_vec())
                                .unwrap()
                        );
                    }
                }
                check_read_content(
                    &buffer_2[header_start..read_end],
                    seq_start - header_start,
                    plus_start - header_start,
                    qual_start - header_start,
                );
            }

            if reporting_level > 0 {
                if run_manager.paired_read_input() {
                    // this is for r1 only if paired end
                    let (qc_total, high_qc) =
                        sum_qc(&buffer_1[qual_start_pr..read_end_pr], check_content);
                    report_manager.update_stats(sample_id, 0, high_qc);
                    report_manager.update_stats(sample_id, 6, qc_total);
                }

                let (qc_total, high_qc) = sum_qc(
                    &buffer_2[read_end - barcode_length..read_end],
                    check_content,
                );
                report_manager.update_stats(sample_id, 2, high_qc);
                report_manager.update_stats(sample_id, 8, qc_total);

                let (qc_total, high_qc) = sum_qc(
                    &buffer_2[qual_start..read_end - barcode_length],
                    check_content,
                );
                report_manager.update_stats(sample_id, shift, high_qc);
                report_manager.update_stats(sample_id, 6 + shift, qc_total);
            }
            //report_manager.update_mismatches(sample_id, 0, 1);
            report_manager.update_mismatches(sample_id, curr_mismatch + 1, 1_u64);

            curr_writing_sample = writing_samples[sample_id];
            // writing preperation
            // this works for mgi format and unde and ambig and ilumina with a bit of extr
            match samples_reads.get_mut(curr_writing_sample) {
                Some(curr_sample) => {
                    curr_sample.compress_and_write(false, &samples_locks[curr_writing_sample]);
                }
                None => {}
            }

            if demultiplex || illumina_format {
                if sample_id >= undetermined_label_id {
                    match samples_reads.get_mut(sample_id) {
                        Some(curr_sample) => {
                            curr_sample.add_barcode_reads(&buffer_2[header_start..read_end + 1]);
                            if run_manager.paired_read_input() {
                                curr_sample
                                    .add_paired_reads(&buffer_1[header_start_pr..read_end_pr + 1]);
                            }
                        }
                        None => {}
                    };
                } else {
                    match samples_reads.get_mut(curr_writing_sample) {
                        Some(curr_sample) => {
                            //curr_sample.add_barcode_reads(&buffer_2[header_start..read_end + 1]);
                            if illumina_format {
                                // Illumina format write the heder and skip the header for mgi.
                                barcode_read_illumina_header_start =
                                    curr_sample.get_barcode_read_compression_end();

                                curr_sample.write_illumina_header(
                                    &buffer_2[header_start..seq_start - raw_shift],
                                    &curr_umi.as_bytes(),
                                    l_position,
                                    sep_position,
                                    !demultiplex,
                                    true,
                                    barcode_read_illumina_header_start == usize::MAX,
                                );
                                if demultiplex {
                                    curr_sample.add_barcode_reads(&curr_barcode.as_bytes());
                                } else if sep_position == seq_start - header_start - 3
                                    && reformated_sample.barcode().len() > 0
                                {
                                    let mut new_barcode = String::from(" ");
                                    new_barcode.push(buffer_2[seq_start - 2] as char);
                                    new_barcode.push_str(":N:0:");
                                    new_barcode.push_str(reformated_sample.barcode());
                                    curr_sample.add_barcode_reads(&new_barcode.as_bytes());
                                }
                                if barcode_read_illumina_header_start < usize::MAX {
                                    curr_sample
                                        .copy_header_2_paired(barcode_read_illumina_header_start);
                                    if run_manager.paired_read_input() {
                                        if demultiplex
                                            || reformated_sample.barcode().len() > 0
                                            || sep_position < seq_start - header_start - 3
                                        {
                                            curr_sample.fix_paired_read_header(
                                                tail_offset,
                                                buffer_1[seq_start_pr - 2 - header_shift],
                                            );
                                        }

                                        header_start_pr = seq_start_pr - 1;
                                    }
                                } else {
                                    if demultiplex {
                                        curr_sample.add_paired_reads(&curr_barcode.as_bytes());
                                    }
                                    if run_manager.paired_read_input() {
                                        curr_sample.fix_paired_read_header(
                                            tail_offset,
                                            buffer_1[seq_start_pr - 2 - header_shift],
                                        );
                                        header_start_pr = seq_start_pr - 1;
                                    }
                                }
                                header_start = seq_start - 1;
                            } else if run_manager.mgi_full_header() {
                                curr_sample
                                    .add_barcode_reads(&buffer_2[header_start..seq_start - 1]);
                                header_info = String::from(" ");
                                if curr_umi.len() > 0 {
                                    header_info.push_str(&curr_umi[1..curr_umi.len()]);
                                    header_info.push(' ');
                                }
                                header_info.push(buffer_2[seq_start - 2] as char);
                                header_info.push_str(&":N:0:");
                                header_info.push_str(&curr_barcode);
                                curr_sample.add_barcode_reads(&header_info.as_bytes());
                                header_start = seq_start - 1;

                                if run_manager.paired_read_input() {
                                    curr_sample.add_paired_reads(
                                        &buffer_1[header_start_pr..seq_start_pr - 1],
                                    );
                                    unsafe {
                                        header_info.as_bytes_mut()[curr_umi.len() + 1] =
                                            buffer_1[seq_start_pr - 2];
                                    }
                                    curr_sample.add_paired_reads(header_info.as_bytes());
                                    header_start_pr = seq_start_pr - 1;
                                }
                            }

                            curr_sample.add_barcode_reads(
                                &buffer_2[header_start..plus_start - writen_barcode_length - 1],
                            );
                            curr_sample.add_barcode_reads(
                                &buffer_2[plus_start - 1..read_end - writen_barcode_length],
                            );
                            curr_sample.add_barcode_reads(&[b'\n']);
                            if run_manager.paired_read_input() {
                                curr_sample
                                    .add_paired_reads(&buffer_1[header_start_pr..read_end_pr + 1]);
                            }
                        }
                        None => {
                            panic!("Read must attached to a specific sample!");
                        }
                    };
                }
            }

            if run_manager.paired_read_input() {
                header_start_pr = read_end_pr + 1;
            }
            header_start = read_end + 1;
            read_cntr += 1;
        }
    }
    return (header_start, header_start_pr);
}

fn analyse_fastq(
    sample_manager: &SampleManager,
    run_manager: &RunManager,
    buffer_info: &BufferInfo,
    all_template_data: Vec<(
        u32,
        HashSet<String>,
        HashSet<String>,
        String,
        HashMap<String, (usize, HashMap<String, usize>)>,
        bool,
        [usize; 10],
    )>,
    reporting_level: usize,
    allowed_mismatches: usize,
    all_index_error: bool,
    buffer_size: usize,
    full_receiver: Receiver<(usize, Vec<u8>, Vec<usize>, usize, Vec<u8>, Vec<usize>)>,
    empty_sender_rb: Sender<(usize, Vec<u8>)>,
    empty_sender_rp: Sender<(usize, Vec<u8>)>,
    parallel_reader: bool,
    samples_locks: Arc<Vec<Mutex<bool>>>,
    demultiplex: bool,
    reformated_sample: &ReformatedSample,
    mut warnings_ls: Vec<bool>,
) -> ReportManager {
    let curr_thread = thread::current().name().unwrap_or("Unnamed").to_string();
    info!("Thread ({}) has started.", curr_thread);
    let mut samples_reads: Vec<SampleData> = if demultiplex {
        create_sample_data_list(
            sample_manager,
            run_manager,
            buffer_info,
            run_manager.read2_has_sequence(),
            run_manager.illumina_format(),
        )
    } else {
        let mut sample = SampleData::new(
            reformated_sample.sample_label().clone(),
            run_manager.paired_read_input(),
            run_manager.read2_has_sequence(),
            buffer_info.clone(),
            run_manager.create_illumina_header_prefix(),
        );
        sample.create_files(
            run_manager.lane(),
            reformated_sample.sample_index(),
            run_manager.illumina_format(),
            run_manager.output_dir(),
            run_manager.paired_read_input(),
        );
        vec![
            sample,
            SampleData::minimal_sample(),
            SampleData::minimal_sample(),
        ]
    };

    //let project_samples = get_project_samples(sample_information).unwrap();
    let mut mismatches_dic_i7: Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>> = Vec::new();
    let mut mismatches_dic_i5: Vec<HashMap<Vec<u8>, (Vec<&String>, usize)>> = Vec::new();

    for template_details in &all_template_data {
        mismatches_dic_i7.push(get_all_mismatches(&template_details.1, allowed_mismatches));
        mismatches_dic_i5.push(get_all_mismatches(&template_details.2, allowed_mismatches));
    }

    let total_samples: usize = sample_manager.get_sample_count();
    let mut report_manager = ReportManager::new(total_samples, allowed_mismatches);
    let paired_input = run_manager.paired_read_input();

    let mut reader_barcode_read = if !parallel_reader {
        Some(get_gzip_reader(run_manager.barcode_reads()))
    } else {
        None
    };

    let mut reader_paired_read = if !parallel_reader && paired_input {
        let r1 = get_gzip_reader(run_manager.paired_reads());
        Some(r1)
    } else {
        None
    };

    let mut main_buffer_1: Vec<u8> = if parallel_reader || !paired_input {
        Vec::new()
    } else {
        vec![0; buffer_size]
    }; // paired read
    let mut main_buffer_2: Vec<u8> = if parallel_reader {
        Vec::new()
    } else {
        vec![0; buffer_size]
    }; // read with barcode

    let mut read_bytes_1: usize = 0;
    let mut read_bytes_2: usize = 0;

    let mut header_start_pr: usize = 0;
    let mut header_start: usize = 0;
    let mut header_start_pr_tmp: usize;
    let mut header_start_tmp: usize;
    //debug!("{}, {}", run_manager.paired_reads().display(), run_manager.barcode_reads().display());
    //debug!("{}, {}  - {} - {}  - {}", buffer_1.len(), buffer_2.len(), minimum_read_bytes, paired_input, parallel_reader);

    //let mut reading_time = Duration::from_secs(0);
    //let start_full = Instant::now();

    loop {
        let (buffer_1, buffer_2, lines_rb, lines_rp) = match parallel_reader {
            true => {
                header_start = 0;
                header_start_pr = 0;
                //let start = Instant::now();
                let (read_bytes_2_tmp, buffer_2, lines_rb, read_bytes_1_tmp, buffer_1, lines_rp) =
                    full_receiver.recv().unwrap();
                read_bytes_2 = read_bytes_2_tmp;
                read_bytes_1 = read_bytes_1_tmp;
                //debug!("Barcode read: Received {} bytes", read_bytes_2);
                //debug!("Paired read: Received {} bytes", read_bytes_1);
                //reading_time += start.elapsed();
                (buffer_1, buffer_2, lines_rb, lines_rp)
            }
            false => {
                //let start = Instant::now();
                read_bytes_2 =
                    read_buffers(read_bytes_2, &mut main_buffer_2, &mut reader_barcode_read);
                if paired_input {
                    read_bytes_1 =
                        read_buffers(read_bytes_1, &mut main_buffer_1, &mut reader_paired_read);
                }
                //reading_time += start.elapsed();
                (Vec::new(), Vec::new(), Vec::new(), Vec::new())
            }
        };

        //debug!("Z1 - Bytes1 {}   -  bytes2   {}", read_bytes_1, read_bytes_2);
        if read_bytes_2 == 0 && read_bytes_1 == 0 {
            info!("Thread ({}) has finished.", curr_thread);
            break;
        } else if (read_bytes_2 == 0 && read_bytes_1 != 0)
            || (read_bytes_2 != 0 && read_bytes_1 == 0 && run_manager.paired_read_input())
        {
            panic!("Something wrong in the input files!");
        }

        (header_start_tmp, header_start_pr_tmp) = process_buffer(
            run_manager,
            buffer_info,
            &mut report_manager,
            if all_index_error {
                allowed_mismatches
            } else {
                allowed_mismatches * 2
            },
            sample_manager,
            &all_template_data,
            &mut samples_reads,
            reporting_level,
            allowed_mismatches,
            if parallel_reader {
                &buffer_1[header_start_pr..read_bytes_1]
            } else {
                &main_buffer_1[header_start_pr..read_bytes_1]
            },
            if parallel_reader {
                &buffer_2[header_start..read_bytes_2]
            } else {
                &main_buffer_2[header_start..read_bytes_2]
            },
            &mismatches_dic_i7,
            &mismatches_dic_i5,
            &samples_locks,
            if parallel_reader {
                lines_rb
            } else {
                memchr_iter(b'\n', &main_buffer_2[header_start..read_bytes_2]).collect()
            },
            if parallel_reader {
                lines_rp
            } else {
                memchr_iter(b'\n', &main_buffer_1[header_start_pr..read_bytes_1]).collect()
            },
            demultiplex,
            reformated_sample,
            &mut warnings_ls,
        );

        //read_leftover_leng_rp
        //debug!("Z1 - 1- {}  -  {}", header_start_pr_tmp, header_start_tmp);
        //debug!("Z1 - 2- {}  -  {}  -  {}  -  {}", read_bytes_1, read_bytes_1 - header_start_pr, read_bytes_2, read_bytes_2 - header_start);
        if parallel_reader {
            //let empty_sender_rb = empty_sender_rb_arc.lock().unwrap();
            match empty_sender_rb.send((0, buffer_2)) {
                Ok(_) => {
                    //debug!("Sending back empty rb after processing");
                }
                Err(e) => println!("2- Error Sending: {}", e),
            }
            if paired_input {
                //let empty_sender_rp = empty_sender_rp_arc.lock().unwrap();
                match empty_sender_rp.send((0, buffer_1)) {
                    Ok(_) => {
                        //debug!("Sending back empty rp after processing");
                    }
                    Err(e) => println!("2- Error Sending: {}", e),
                }
            }
        } else {
            header_start += header_start_tmp;
            header_start_pr += header_start_pr_tmp;
            if read_bytes_2 - header_start < 10000 {
                //info!("ZB0: {} - {} - *{}*", header_start, read_bytes_2 - header_start, String::from_utf8(buffer_2[header_start..read_bytes_2].to_vec()).unwrap());
                copy_within_a_slice(
                    &mut main_buffer_2,
                    header_start,
                    0,
                    read_bytes_2 - header_start,
                );
                read_bytes_2 -= header_start;
                header_start = 0;
                //info!("ZB00: {} - {} - *{}*", header_start, read_bytes_2 - header_start, String::from_utf8(buffer_2[header_start..read_bytes_2].to_vec()).unwrap());
                //prv_h_s = header_start;
                //prv_end = read_bytes_2;
            }

            if paired_input {
                if read_bytes_1 - header_start_pr < 10000 {
                    copy_within_a_slice(
                        &mut main_buffer_1,
                        header_start_pr,
                        0,
                        read_bytes_1 - header_start_pr,
                    );
                    read_bytes_1 -= header_start_pr;
                    header_start_pr = 0;
                }
            }
        }
        //debug!("Z1 - 3- {}  -  {}  -  {}  -  {}", read_bytes_1, header_start_pr, read_bytes_2, header_start);
    }

    //debug!("Spent {:?} for reading", reading_time);
    //debug!("Spent {:?} for evertyihng else", start_full.elapsed() - reading_time);

    for curr_writing_sample in 0..total_samples {
        if let Some(curr_sample) = samples_reads.get_mut(curr_writing_sample) {
            curr_sample.compress_and_write(true, &samples_locks[curr_writing_sample]);
        };
    }

    // should update run reports on here
    report_manager.set_sample_total_reads();
    report_manager
}

pub fn initiate_demultiplexing(demultiplex_command: &ArgMatches) {
    //let arg_input_folder_path: &String = demultiplex_command.get_one::<String>("arg_input_folder_path").unwrap();
    //let mut arg_read1_file_path: String = demultiplex_command.get_one::<String>("arg_read1_file_path").unwrap().to_string();
    //let mut arg_read2_file_path: String = demultiplex_command.get_one::<String>("arg_read2_file_path").unwrap().to_string();
    //let arg_ouput_dir: &String = demultiplex_command.get_one::<String>("arg_ouput_dir").unwrap();
    //let arg_report_dir: &String = demultiplex_command.get_one::<String>("arg_report_dir").unwrap();
    //let arg_lane: &String = demultiplex_command.get_one::<String>("arg_lane").unwrap();
    //let arg_instrument: &String = demultiplex_command.get_one::<String>("arg_instrument").unwrap();
    //let arg_run: &String = demultiplex_command.get_one::<String>("arg_run").unwrap();
    //let arg_force: &bool = demultiplex_command.get_one::<bool>("arg_force").unwrap();
    //let arg_not_mgi: &bool = demultiplex_command.get_one::<bool>("arg_not_mgi").unwrap();
    //let arg_read1_file_name_suf: &String = demultiplex_command.get_one::<String>("arg_read1_file_name_suf").unwrap();
    //let arg_read2_file_name_suf: &String = demultiplex_command.get_one::<String>("arg_read2_file_name_suf").unwrap();
    //let arg_info_file: &String = demultiplex_command.get_one::<String>("arg_info_file").unwrap();
    let sample_manager = SampleManager::new(
        demultiplex_command
            .get_one::<String>("arg_sample_sheet_file_path")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_template")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<bool>("arg_i7_rc")
            .unwrap()
            .clone(),
        demultiplex_command
            .get_one::<bool>("arg_i5_rc")
            .unwrap()
            .clone(),
        demultiplex_command
            .get_one::<String>("arg_undetermined_label")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_ambiguous_label")
            .unwrap()
            .to_string(),
    );

    let mut run_manager = RunManager::new(
        demultiplex_command
            .get_one::<String>("arg_input_folder_path")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_read1_file_path")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_read2_file_path")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_ouput_dir")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_report_dir")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_lane")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_instrument")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_run")
            .unwrap()
            .to_string(),
        !*demultiplex_command.get_one::<bool>("arg_not_mgi").unwrap(),
        *demultiplex_command.get_one::<bool>("arg_force").unwrap(),
        demultiplex_command
            .get_one::<String>("arg_read1_file_name_suf")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_read2_file_name_suf")
            .unwrap()
            .to_string(),
        demultiplex_command
            .get_one::<String>("arg_info_file")
            .unwrap()
            .to_string(),
        !*demultiplex_command
            .get_one::<bool>("arg_disable_illumina_format")
            .unwrap(),
        *demultiplex_command
            .get_one::<bool>("arg_keep_barcode")
            .unwrap(),
        *demultiplex_command
            .get_one::<bool>("arg_comprehensive_scan")
            .unwrap(),
        *demultiplex_command
            .get_one::<bool>("arg_ignore_undetermined")
            .unwrap(),
        *demultiplex_command
            .get_one::<bool>("arg_check_content")
            .unwrap(),
        *demultiplex_command
            .get_one::<bool>("arg_mgi_full_header")
            .unwrap(),
    );
    let (barcode_read_info, paired_read_info) = run_manager.get_read_information();
    if run_manager.lane().len() == 0 {
        info!("lane detected in the read header will be used for this run!");
        run_manager.set_lane(barcode_read_info.lane().clone());
    }
    let barcode_length: usize = sample_manager.all_template_data()[0].6[9];
    debug!(
        "barcode length = {}, read2 length = {}",
        barcode_length,
        barcode_read_info.sequence_length()
    );
    if barcode_length == *barcode_read_info.sequence_length() {
        info!("It is assumed that read 2 contains barcode only without read sequence!");
        run_manager.set_read2_has_sequence(false);
    }

    run_manager.confirm_format();
    let compression_buffer_size = *demultiplex_command
        .get_one::<usize>("arg_compression_buffer_size")
        .unwrap();
    let writing_buffer_size = *demultiplex_command
        .get_one::<usize>("arg_writing_buffer_size")
        .unwrap();

    let arg_allowed_mismatches: &usize = demultiplex_command
        .get_one::<usize>("arg_allowed_mismatches")
        .unwrap();
    let mut buffer_info = BufferInfo::new(
        writing_buffer_size,
        compression_buffer_size,
        *demultiplex_command
            .get_one::<u32>("arg_compression_level")
            .unwrap(),
        compression_buffer_size - barcode_read_info.read_length() - paired_read_info.read_length(),
        writing_buffer_size,
    );
    let arg_memory: f64 = *demultiplex_command.get_one::<f64>("arg_memory").unwrap();
    if 0.0 < arg_memory && arg_memory <= 0.5 {
        panic!("Requested memory should be greater than 0.5 GB!");
    }

    let mut tmp_reader_threads = *demultiplex_command
        .get_one::<usize>("arg_threads_r")
        .unwrap();

    let mut tmp_processing_threads = *demultiplex_command
        .get_one::<usize>("arg_threads_w")
        .unwrap();

    let (reader_threads, processing_threads) = if tmp_processing_threads > 0
        && tmp_reader_threads > 0
    {
        if run_manager.paired_read_input() {
            if tmp_reader_threads > 4 {
                warn!("Reader threads can not be more than 4! extra threads will be used for processing!");
                tmp_processing_threads += tmp_reader_threads - 4;
                tmp_reader_threads = 4;
            } else if tmp_reader_threads == 3 {
                warn!("Reader threads shoyuld be either 0, 1, 2, or 4! extra threads will be used for processing!");
                tmp_reader_threads = 2;
                tmp_processing_threads += 1;
            }
        } else {
            if tmp_reader_threads > 2 {
                warn!("Reader threads can not be more than 2 for single end input!, extra threads will be used for processing!");
                tmp_processing_threads += tmp_reader_threads - 2;
                tmp_reader_threads = 2;
            }
        }
        info!(
            "Reader threads ({}) and processing threads ({}) will be used!",
            tmp_reader_threads, tmp_processing_threads
        );
        (tmp_reader_threads, tmp_processing_threads)
    } else {
        get_cpus(
            *demultiplex_command.get_one::<usize>("arg_threads").unwrap(),
            run_manager.paired_read_input(),
        )
    };

    let max_buffer_size = calculate_largest_buffer_size(
        get_available_memory(arg_memory),
        sample_manager.get_sample_count(),
        buffer_info.compression_buffer_size(),
        !run_manager.paired_read_input(),
        processing_threads,
    );
    buffer_info.calculate_final_writing_buffer_size(max_buffer_size);

    //let arg_writing_buffer_size: &usize = demultiplex_command.get_one::<usize>("arg_writing_buffer_size").unwrap();
    let arg_report_limit: &usize = demultiplex_command
        .get_one::<usize>("arg_report_limit")
        .unwrap();
    let arg_report_level: &usize = demultiplex_command
        .get_one::<usize>("arg_report_level")
        .unwrap();
    info!("Reporting level is: {}", arg_report_level);
    //let arg_compression_level: &u32 = demultiplex_command.get_one::<u32>("arg_compression_level").unwrap();
    //let arg_compression_buffer_size:  &usize = demultiplex_command.get_one::<usize>("arg_compression_buffer_size").unwrap();
    let arg_all_index_error: &bool = demultiplex_command
        .get_one::<bool>("arg_all_index_error")
        .unwrap();

    match demultiplex(
        &sample_manager,
        &run_manager,
        &buffer_info,
        *arg_allowed_mismatches,
        *arg_report_level,
        *arg_all_index_error,
        reader_threads,
        processing_threads,
    ) {
        Ok(mut report_manager) => {
            let max_mismatches = if *arg_all_index_error {
                *arg_allowed_mismatches + 1
            } else {
                *arg_allowed_mismatches * 2 + 1
            };
            let shift = if run_manager.paired_read_input() {
                1
            } else {
                0
            };
            report_manager.prepare_final_data(
                max_mismatches,
                shift,
                &barcode_read_info,
                &paired_read_info,
                barcode_length,
            );
            report_manager.write_reports(
                &run_manager,
                &sample_manager,
                *arg_report_level,
                *arg_report_limit,
                max_mismatches,
                usize::MAX,
            );
        }
        Err(err) => eprintln!("Error: {}", err),
    };
}

pub fn merge_qc_reports(
    qc_report_paths: &[String],
    output_dir: &String,
    lane: &String,
    project: &String,
) {
    if qc_report_paths.len() == 0 {
        panic!("report directories are not provided!");
    }

    if output_dir.len() == 0 {
        info!("output will be written to the current work directory!");
    }
    let used_lane = if lane.len() == 0 {
        String::from("all")
    } else {
        lane.to_string()
    };

    let mut project_id: String;
    let mut sample_id: String;
    let mut vals: Vec<String>;
    let mut report_manager = ReportManager::new(1000, 10);
    let mut sample_manager = SampleManager::default();
    let flowcell_id = if project.len() == 0 {
        Path::new(&qc_report_paths[0])
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .split(".")
            .collect::<Vec<&str>>()[0]
            .to_string()
    } else {
        project.to_string()
    };
    let mut run_manager = RunManager::default();
    run_manager.set_flowcell(flowcell_id);
    run_manager.set_lane(used_lane);
    run_manager.set_report_dir(PathBuf::from(output_dir));

    let mut sample_index;
    let mut max_mismatches = 0;

    for qc_report_path in qc_report_paths {
        info!("Reading {} ...", qc_report_path);
        let file_content = fs::read_to_string(Path::new(qc_report_path)).unwrap();
        let lines = file_content.lines();

        for line in lines {
            if line.starts_with("job_number\tsample_id") {
                continue;
            }
            //println!("{}", line);
            vals = line.split("\t").map(|x| x.to_string()).collect();
            project_id = vals.remove(0);
            sample_id = vals.remove(0);
            sample_index = sample_manager.get_sample_index(&sample_id);
            if sample_index == usize::MAX {
                sample_manager.add_sample(&sample_id, &project_id);
                sample_index = sample_manager.get_sample_count();
                report_manager.add_stats_entry();
                report_manager.add_mismatches_entry(20);
            }
            sample_manager.add_project_sample(&project_id, sample_index);
            if vals.len() - 10 > max_mismatches {
                max_mismatches = vals.len() - 10;
            }
            for i in 0..10 {
                report_manager.update_stats(sample_index, i, vals[i].parse::<u64>().unwrap());
            }
            for i in 9..vals.len() {
                report_manager.update_mismatches(sample_index, i, vals[i].parse::<u64>().unwrap());
            }
        }
    }
    report_manager.truncate_mismatches(max_mismatches);

    report_manager.write_reports(
        &run_manager,
        &sample_manager,
        2,
        0,
        max_mismatches,
        usize::MAX,
    );

    //for (sample_id, val) in map.iter_mut() {  }
}

pub fn detect_template(template_command: &ArgMatches) {
    let start = Instant::now();
    let sample_manager = SampleManager::from_simple_sheet(
        template_command
            .get_one::<String>("arg_sample_sheet_file_path")
            .unwrap()
            .to_string(),
    );

    let path = Path::new(template_command.get_one::<String>("arg_ouput_dir").unwrap());

    let output_file = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_string();

    let parent_path = path.parent().and_then(|p| p.to_str()).unwrap_or("");

    let mut run_manager = RunManager::new(
        template_command
            .get_one::<String>("arg_input_folder_path")
            .unwrap()
            .to_string(),
        template_command
            .get_one::<String>("arg_read1_file_path")
            .unwrap()
            .to_string(),
        template_command
            .get_one::<String>("arg_read2_file_path")
            .unwrap()
            .to_string(),
        parent_path.to_string(),
        String::new(),
        String::new(),
        String::new(),
        String::new(),
        true,
        true,
        template_command
            .get_one::<String>("arg_read1_file_name_suf")
            .unwrap()
            .to_string(),
        template_command
            .get_one::<String>("arg_read2_file_name_suf")
            .unwrap()
            .to_string(),
        template_command
            .get_one::<String>("arg_info_file")
            .unwrap()
            .to_string(),
        true,
        false,
        true,
        true,
        true,
        false,
    );

    let input_barcode_length: usize = *template_command
        .get_one::<usize>("arg_barcode_length")
        .unwrap();
    let testing_reads: usize = *template_command
        .get_one::<usize>("arg_testing_reads")
        .unwrap();
    let max_umi_length: usize = *template_command
        .get_one::<usize>("arg_max_umi_length")
        .unwrap();
    let use_popular_template: bool = *template_command
        .get_one::<bool>("arg_popular_template")
        .unwrap();
    let add_umi: bool = !template_command.get_one::<bool>("arg_no_umi").unwrap();

    let (barcode_read_info, paired_read_info) = run_manager.get_read_information();
    let sample_indexes = sample_manager.get_samples_indices().unwrap();

    let barcode_length: usize;
    if input_barcode_length > 0 {
        barcode_length = input_barcode_length as usize;
    } else {
        barcode_length = barcode_read_info.sequence_length() - paired_read_info.sequence_length();
        info!("Barcode length is calculated as the difference between R2 length and R1.");
        let max_barcode_length: usize = match sample_indexes
            .iter()
            .map(|it| it[0].len() + it[2].len())
            .max()
        {
            Some(tmp_max) => tmp_max,
            None => panic!("Sample sheet should have samples!"),
        };

        if barcode_length > max_barcode_length + max_umi_length {
            panic!(
                "The difference in read length is {}. It is greater than the the length of indexes and possible UMI {}. You need to prvide barcode length for this run!",
                barcode_length,
                max_barcode_length + max_umi_length
            );
        }
    }

    info!("Barcode length: {}", barcode_length);

    let mut read_barcode_seq;
    let mut read_bytes;

    let mut read_cntr: usize = 0;
    let mut matches_stat: HashMap<String, Vec<usize>> = HashMap::new();
    let sample_information = sample_manager.sample_information();
    let mut barcode_reader = get_buf_reader(run_manager.barcode_reads());

    loop {
        if read_cntr >= testing_reads {
            break;
        }
        read_barcode_seq = String::new();
        read_bytes = barcode_reader.read_line(&mut read_barcode_seq).unwrap();
        if read_bytes == 0 {
            break;
        }
        read_barcode_seq = String::new();
        barcode_reader.read_line(&mut read_barcode_seq).unwrap();
        read_barcode_seq = read_barcode_seq
            [read_barcode_seq.len() - barcode_length - 1..read_barcode_seq.len() - 1]
            .to_string();
        //println!(" {} Barcode seq: {}", read_cntr, read_barcode_seq);

        for sample_itr in 0..sample_indexes.len() {
            let sample_i7 = &sample_indexes[sample_itr][0];
            let sample_i7_rc = &sample_indexes[sample_itr][1];
            let sample_i5 = &sample_indexes[sample_itr][2];
            let sample_i5_rc = &sample_indexes[sample_itr][3];

            find_matches(
                sample_itr,
                &mut matches_stat,
                sample_indexes.len(),
                &read_barcode_seq,
                &sample_i7,
                &sample_i7_rc,
                &sample_i5,
                &sample_i5_rc,
            );
        }

        barcode_reader.read_line(&mut read_barcode_seq).unwrap();
        barcode_reader.read_line(&mut read_barcode_seq).unwrap();
        read_cntr += 1;
    }

    let mut sample_reads_final: Vec<Vec<(String, usize)>> = Vec::new();

    for sample_itr in 0..sample_information.len() {
        let mut no_matches = true;
        sample_reads_final.push(Vec::new());
        for (template_str, sample_reads) in &matches_stat {
            if sample_reads[sample_itr] > 0 {
                sample_reads_final[sample_itr]
                    .push((template_str.to_string(), sample_reads[sample_itr]));
                no_matches = false;
            }
        }
        if no_matches {
            warn!(
                "Sample ({}) has no matches with any possible template! might be an issue with its indexes!",
                sample_information[sample_itr][SAMPLE_COLUMN].to_string()
            );
            warn!("Consder increasing the '--testing-reads' parameter if the indexes are correct!");
        }
        sample_reads_final[sample_itr].sort_by(|a, b| b.1.cmp(&a.1));
    }

    let mut popular_template = String::new();
    let mut popular_template_cnt = 0;

    for (template_str, sample_reads) in &matches_stat {
        let cnt = sample_reads.iter().sum();
        if cnt > popular_template_cnt {
            popular_template = template_str.to_string();
            popular_template_cnt = cnt;
        }
    }

    if popular_template_cnt == 0 {
        panic!("Something wrong, there is no matches with any sample!");
    }

    if use_popular_template {
        let tmp = get_mgikit_template(
            &popular_template,
            add_umi,
            barcode_length,
            sample_indexes[0][0].len(),
            sample_indexes[0][2].len(),
        );
        info!(
            "The most frequent template is {} - Appeared {} times!",
            format!("{} i7_rc is {} and i5_rc is {}", tmp.0, tmp.1, tmp.2),
            popular_template_cnt
        );
    }

    let mut out_str = String::from("sample_id\ti7\ti5\tjob_number\ttemplate\ti7_rc\ti5_rc\n");

    let mut out_str_full = String::from("sample_id\ti7\ti5\tall-matches\n");

    for sample_itr in 0..sample_information.len() {
        out_str.push_str(&sample_information[sample_itr][SAMPLE_COLUMN]);
        out_str.push('\t');

        out_str.push_str(&sample_information[sample_itr][I7_COLUMN]);
        out_str.push('\t');

        out_str.push_str(&sample_information[sample_itr][I5_COLUMN]);
        out_str.push('\t');

        out_str.push_str(&sample_information[sample_itr][PROJECT_ID_COLUMN]);
        out_str.push('\t');

        if use_popular_template {
            if sample_reads_final[sample_itr].len() > 0
                && popular_template != sample_reads_final[sample_itr][0].0
            {
                info!(
                    "Using popular template for sample {} instead of its detected template {}!",
                    sample_information[sample_itr][SAMPLE_COLUMN],
                    sample_reads_final[sample_itr][0].0
                );
            } else if sample_reads_final[sample_itr].len() == 0 {
                info!(
                    "Using popular template for sample {} as no matches were found!",
                    sample_information[sample_itr][SAMPLE_COLUMN]
                );
            }

            let template_info = get_mgikit_template(
                &popular_template,
                add_umi,
                barcode_length,
                sample_information[sample_itr][I7_COLUMN].len(),
                sample_information[sample_itr][I5_COLUMN].len(),
            );
            out_str.push_str(&template_info.0);
            out_str.push('\t');
            out_str.push_str(&template_info.1);
            out_str.push('\t');
            out_str.push_str(&template_info.2);
        } else {
            let template_info = match sample_reads_final[sample_itr].len() > 0 {
                true => get_mgikit_template(
                    &sample_reads_final[sample_itr][0].0,
                    add_umi,
                    barcode_length,
                    sample_indexes[sample_itr][0].len(),
                    sample_indexes[sample_itr][2].len(),
                ),
                false => (
                    String::from("No matches with this sample"),
                    String::from("."),
                    String::from("."),
                ),
            };
            out_str.push_str(&template_info.0);
            out_str.push('\t');
            out_str.push_str(&template_info.1);
            out_str.push('\t');
            out_str.push_str(&template_info.2);
        }

        out_str.push('\n');

        out_str_full.push_str(&sample_information[sample_itr][SAMPLE_COLUMN]);
        out_str_full.push('\t');

        out_str_full.push_str(&sample_information[sample_itr][I7_COLUMN]);
        out_str_full.push('\t');

        out_str_full.push_str(&sample_information[sample_itr][I5_COLUMN]);
        out_str_full.push('\t');

        if sample_reads_final[sample_itr].len() > 0 {
            let tmp: Vec<String> = sample_reads_final[sample_itr]
                .iter()
                .map(|tmp_det| format!("{} ({})", tmp_det.0, tmp_det.1))
                .collect();
            out_str_full.push_str(&tmp.join("\t"));
        } else {
            out_str_full.push_str(&String::from("No matches found!"));
        }

        out_str_full.push('\n');
    }

    write_file(
        &format!(
            "{}/{}_template.tsv",
            run_manager.output_dir().display(),
            output_file
        ),
        &out_str,
    );
    write_file(
        &format!(
            "{}/{}_details.tsv",
            run_manager.output_dir().display(),
            output_file
        ),
        &out_str_full,
    );

    let dur = start.elapsed();

    info!(
        "{} reads were processed in {} secs.",
        read_cntr,
        dur.as_secs()
    );
}

pub fn reformat(reformat_command: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let reporting_level: usize = *reformat_command
        .get_one::<usize>("arg_report_level")
        .unwrap();

    let start = Instant::now();
    let dur: std::time::Duration;

    let mut run_manager = RunManager::new(
        String::new(),
        reformat_command
            .get_one::<String>("arg_read1_file_path")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_read2_file_path")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_ouput_dir")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_report_dir")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_lane")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_instrument")
            .unwrap()
            .to_string(),
        reformat_command
            .get_one::<String>("arg_run")
            .unwrap()
            .to_string(),
        true,
        *reformat_command.get_one::<bool>("arg_force").unwrap(),
        String::new(),
        String::new(),
        reformat_command
            .get_one::<String>("arg_info_file")
            .unwrap()
            .to_string(),
        !*reformat_command
            .get_one::<bool>("arg_disable_illumina_format")
            .unwrap(),
        false,
        false,
        false,
        *reformat_command
            .get_one::<bool>("arg_check_content")
            .unwrap(),
        false,
    );
    let (barcode_read_info, paired_read_info) = run_manager.get_read_information();
    if run_manager.lane().len() == 0 {
        info!("lane detected in the read header will be used for this run!");
        run_manager.set_lane(barcode_read_info.lane().clone());
    }

    run_manager.confirm_format();
    let compression_buffer_size = *reformat_command
        .get_one::<usize>("arg_compression_buffer_size")
        .unwrap();
    let writing_buffer_size = *reformat_command
        .get_one::<usize>("arg_writing_buffer_size")
        .unwrap();

    let (sample_label_sb, _, sb_lane, _) = parse_sb_file_name(
        &run_manager
            .barcode_reads()
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string(),
    );

    let mut sample_label = reformat_command
        .get_one::<String>("arg_sample_label")
        .unwrap()
        .to_string();
    if sample_label.len() == 0 {
        if sample_label_sb.len() == 0 {
            panic!(
                "Sample label needs to be passed either through '--sample-label' parameter or the third part of the file name after separation by '_'"
            );
        }
        sample_label = sample_label_sb;
    }
    if sb_lane != *run_manager.lane() {
        warn!(
            "The lane extracted from the file name ({}) is not the same as the provided lane ({})",
            sb_lane,
            run_manager.lane()
        );
    }
    let reformated_sample = ReformatedSample::new(
        sample_label.clone(),
        *reformat_command
            .get_one::<usize>("arg_sample_index")
            .unwrap(),
        *reformat_command.get_one::<usize>("arg_umi_length").unwrap(),
        reformat_command
            .get_one::<String>("arg_barcode")
            .unwrap()
            .to_string(),
    );

    let mut buffer_info = BufferInfo::new(
        writing_buffer_size,
        compression_buffer_size,
        *reformat_command
            .get_one::<u32>("arg_compression_level")
            .unwrap(),
        compression_buffer_size - barcode_read_info.read_length() - paired_read_info.read_length(),
        writing_buffer_size,
    );

    let arg_memory: f64 = *reformat_command.get_one::<f64>("arg_memory").unwrap();
    if 0.0 < arg_memory && arg_memory <= 0.5 {
        panic!("Requested memory should be greater than 0.5 GB!");
    }
    let max_buffer_size = calculate_largest_buffer_size(
        get_available_memory(arg_memory),
        1,
        buffer_info.compression_buffer_size(),
        !run_manager.paired_read_input(),
        1,
    );
    buffer_info.calculate_final_writing_buffer_size(max_buffer_size);

    if run_manager.force() {
        let mut sample = SampleData::new(
            sample_label.clone(),
            run_manager.paired_read_input(),
            run_manager.read2_has_sequence(),
            buffer_info.clone(),
            run_manager.create_illumina_header_prefix(),
        );
        sample.create_files(
            run_manager.lane(),
            reformated_sample.sample_index(),
            run_manager.illumina_format(),
            run_manager.output_dir(),
            run_manager.paired_read_input(),
        );
        sample.delete_sample_files();
    }

    let sample_lockes = vec![Mutex::new(false), Mutex::new(false), Mutex::new(false)];
    let samples_locks_arc = Arc::new(sample_lockes);
    let (empty_sender_dummy, _) = bounded(1);
    let (_, full_receiver_dummy) = bounded(1);
    let sample_manager = SampleManager::dummy_sample(sample_label.clone());
    info!("Reporting level is: {}", reporting_level);
    let mut report_manager = analyse_fastq(
        &sample_manager,
        &run_manager,
        &buffer_info,
        Vec::new(),
        reporting_level,
        5,
        false,
        BUFFER_SIZE,
        full_receiver_dummy,
        empty_sender_dummy.clone(),
        empty_sender_dummy,
        false,
        samples_locks_arc,
        false,
        &reformated_sample,
        vec![true; 10],
    );

    if reporting_level > 0 {
        report_manager.prepare_final_data(
            5,
            if run_manager.paired_read_input() {
                1
            } else {
                0
            },
            &barcode_read_info,
            &paired_read_info,
            reformated_sample.umi_length(),
        );
    }
    report_manager.write_reports(&run_manager, &sample_manager, reporting_level, 0, 5, 0);
    dur = start.elapsed();
    info!(
        "{} reads were processed in {} secs.",
        report_manager.get_total_reads(),
        dur.as_secs()
    );
    Ok(())
}
