use getset::{ Getters, Setters, CopyGetters };
use crate::file_utils::*;
use chrono::prelude::Local;
use std::fs::File;
use std::io::{ self, BufRead };
use log::{ info, error, warn };
use core::panic;
use std::fs;
use std::path::Path;
use std::path::PathBuf;

#[derive(Getters, Setters, CopyGetters, Clone, Default)]
pub struct RunManager {
    #[getset(get = "pub", set = "pub")]
    input_dir: PathBuf,
    #[getset(get = "pub", set = "pub")]
    barcode_reads: PathBuf,
    #[getset(get = "pub", set = "pub")]
    paired_reads: PathBuf,
    #[getset(get = "pub", set = "pub")]
    output_dir: PathBuf,
    #[getset(get = "pub", set = "pub")]
    report_dir: PathBuf,
    #[getset(get = "pub", set = "pub")]
    lane: String,
    #[getset(get = "pub", set = "pub")]
    instrument: String,
    #[getset(get = "pub", set = "pub")]
    run: String,
    #[getset(get_copy = "pub", set = "pub")]
    mgi_data: bool,
    #[getset(get_copy = "pub", set = "pub")]
    force: bool,
    #[getset(get_copy = "pub", set = "pub")]
    illumina_format: bool,
    #[getset(get_copy = "pub", set = "pub")]
    keep_barcode: bool,
    #[getset(get_copy = "pub", set = "pub")]
    comprehensive_scan: bool,
    #[getset(get_copy = "pub", set = "pub")]
    ignore_undetermined: bool,
    #[getset(get = "pub", set = "pub")]
    barcode_read_info: ReadInfo,
    #[getset(get = "pub", set = "pub")]
    paired_read_info: ReadInfo,
    #[getset(get_copy = "pub", set = "pub")]
    read2_has_sequence: bool,
    #[getset(get_copy = "pub", set = "pub")]
    check_content: bool,
    #[getset(get_copy = "pub", set = "pub")]
    mgi_full_header: bool,
}

impl RunManager {
    pub fn new(
        input_dir: String,
        read1_file_path: String,
        read2_file_path: String,
        ouput_dir: String,
        report_dir: String,
        mut lane: String,
        instrument: String,
        run: String,
        mgi_data: bool,
        force: bool,
        read1_file_name_suf: String,
        read2_file_name_suf: String,
        info_file: String,
        illumina_format: bool,
        keep_barcode: bool,
        comprehensive_scan: bool,
        ignore_undetermined: bool,
        check_content: bool,
        mgi_full_header: bool
    ) -> Self {
        let (paired_read_file_path_final, read_barcode_file_path_final, _) = if input_dir.len() > 0 {
            info!("Input directory: {}", &input_dir);
            get_read_files_from_input_dir(&input_dir, &read1_file_name_suf, &read2_file_name_suf)
        } else {
            validate_and_assigne_input_reads(&read1_file_path, &read2_file_path)
        };

        let (mut final_instrument, mut final_run) = if instrument.len() == 0 || run.len() == 0 {
            parse_info_file(
                &find_info_file(
                    &PathBuf::from(info_file),
                    &PathBuf::new(),
                    &read_barcode_file_path_final
                )
            )
        } else {
            (String::new(), String::new())
        };

        if instrument.len() > 0 {
            final_instrument = instrument.clone();
        }
        if run.len() > 0 {
            final_run = run.clone();
        }

        if lane.len() == 0 {
            lane = get_lane_from_file(&read_barcode_file_path_final, 1, '_');
        }

        let (output_directory, report_directory) = prepare_output_report_dir(
            &ouput_dir,
            &report_dir,
            force
        );
        info!("Output directory: {}", output_directory.display());
        info!("Reports directory: {}", report_directory.display());
        info!("MGI input fastq files:  {}", mgi_data);
        info!("Validate fastq content: {}", check_content);
        info!("Trim Barcode: {}", !keep_barcode);
        if illumina_format && !mgi_data {
            panic!(
                "mgikit does not refomat output files in Illumina foramt unless the input fastq files are in MGI format! Disable `--not-mgi` or enable `--disable-illumina-format`"
            );
        }

        Self {
            input_dir: PathBuf::from(input_dir),
            barcode_reads: read_barcode_file_path_final,
            paired_reads: paired_read_file_path_final,
            output_dir: output_directory,
            report_dir: report_directory,
            lane,
            instrument: final_instrument,
            run: final_run,
            mgi_data,
            force,
            illumina_format,
            keep_barcode,
            comprehensive_scan,
            ignore_undetermined,
            read2_has_sequence: true,
            barcode_read_info: ReadInfo::default(),
            paired_read_info: ReadInfo::default(),
            check_content,
            mgi_full_header,
        }
    }

    pub fn confirm_format(&self) {
        if self.illumina_format {
            info!("Output format is Illumina");
            if self.run.is_empty() || self.instrument.is_empty() || self.flowcell().is_empty() {
                error!(
                    "Missing information for Illumina header! instrument: {}, lane: {}, flowcell: {}",
                    self.instrument,
                    self.lane,
                    self.flowcell()
                );
                panic!("Missing mandatory information!");
            }
        } else {
            info!("Output format is MGI");
            if self.mgi_full_header {
                info!("Sample barcode and UMI if available will be written into the read header.");
            }
        }
    }

    pub fn create_illumina_header_prefix(&self) -> String {
        if self.instrument.len() == 0 || self.lane.len() == 0 || self.flowcell().len() == 0 {
            error!(
                "Missing information for Illumina header! instrument: {}, lane: {}, flowcell: {}",
                self.instrument,
                self.lane,
                self.flowcell()
            );
            panic!("Missing mandatory information!");
            //return String::new();
        }
        let mut header = String::from("@");
        header.push_str(&self.instrument);
        header.push(':');
        header.push_str(&self.run);
        header.push(':');
        header.push_str(&self.flowcell());
        header.push(':');
        header
    }

    pub fn l_position(&self) -> usize {
        self.barcode_read_info.l_position().clone()
    }

    pub fn paired_read_input(&self) -> bool {
        self.paired_read_info.read_length > 0
    }

    pub fn flowcell(&self) -> String {
        self.barcode_read_info.flowcell.to_string()
    }

    pub fn barcode_read_sequence_len(&self) -> usize {
        self.barcode_read_info.sequence_length().clone()
    }

    pub fn paired_read_sequence_len(&self) -> usize {
        self.paired_read_info.sequence_length().clone()
    }

    pub fn barcode_read_len(&self) -> usize {
        self.barcode_read_info.read_length().clone()
    }

    pub fn paired_read_len(&self) -> usize {
        self.paired_read_info.read_length().clone()
    }

    pub fn set_flowcell(&mut self, flowcell: String) {
        self.barcode_read_info.set_flowcell(flowcell);
    }

    pub fn get_read_information(&mut self) -> (ReadInfo, ReadInfo) {
        /*
        if self.lane.len() == 0 {
            lane = format!("L0{}", header_lane);
        }else{
            if lane != format!("L0{}", header_lane){
                warn!("The lane in the read header (L0{}) does not match with the lane provided or extracted from the input file name {}!", header_lane, lane);
            }
        }
        if illumina_format {
            info!("Read header and Output files: Illumina format.");
            info!("Instrumnet: {}", instrument);
            info!("Run: {}", run);
            info!("Lane: {}", lane);
            if lane.len() == 0 || instrument.len() == 0 || run.len() == 0 {
                panic!("Instrument id, lane number and run number are required for QC reports and when Illumina format is requested!")
            }
        }else {
            info!("Read header and Output files: MGI format.");
            info!("Lane: {}", lane);
            if lane.len() == 0 {
                panic!("Lane number is required for QC reports!")
            }
        }
        */

        //  Get reads information
        let whole_paired_read_len: usize;
        let barcode_read_length: usize;
        let mut paired_read_length: usize = 0;
        let only_plus_r1: bool;
        let flowcell: String;
        let header_lane: String;
        let l_position: usize;
        let only_plus_r2: bool;
        let dynamic_demultiplexing = false;

        let mut reader_barcode_read_tmp = get_buf_reader(&self.barcode_reads);
        let (header, seq, plus, quality) = get_read_parts(&mut reader_barcode_read_tmp);

        let whole_read_barcode_len = header.len() + seq.len() + plus.len() + quality.len();

        if !self.mgi_data {
            flowcell = Local::now().format("%Y%m%dT%H%M%S").to_string();
            l_position = 0;
            header_lane = String::new();
        } else {
            (flowcell, l_position, header_lane) = get_flowcell_lane_info(&header);
            //info!("Detected flowcell from the header of the first read is {}.", flowcell);
            info!("Detected lane from the header of the first read is {}.", header_lane);
            if !"1234".contains(&header_lane) {
                warn!("The detected lane ( = {}) is not recognised! Expected 1, 2, 3 or 4!", header_lane);
            }
        }
        barcode_read_length = seq.chars().count() - 1;
        only_plus_r2 = plus == "+\n";
        let barcode_read = ReadInfo::new(
            barcode_read_length,
            whole_read_barcode_len,
            l_position,
            flowcell,
            format!("L0{}", header_lane)
        );

        if !only_plus_r2 && !dynamic_demultiplexing {
            panic!(
                "Expected read format is not satisified. You can try rerunning using --flexible parameter."
            );
        }

        let paired_read = if self.paired_reads().exists() {
            let mut reader_paired_read_buff = get_buf_reader(&self.paired_reads);
            let (header, seq, plus, quality) = get_read_parts(&mut reader_paired_read_buff);
            whole_paired_read_len = header.len() + seq.len() + plus.len() + quality.len();
            //header_length_r1 = header.len();
            paired_read_length = seq.chars().count() - 1;
            only_plus_r1 = plus == "+\n";
            if !only_plus_r1 {
                panic!(
                    "Expected read format is not satisified. You can try running --flexible command."
                );
            }
            ReadInfo::new(
                paired_read_length,
                whole_paired_read_len,
                0,
                String::new(),
                String::new()
            )
        } else {
            ReadInfo::new(0, 0, 0, String::new(), String::new())
        };
        self.paired_read_info = paired_read.clone();
        self.barcode_read_info = barcode_read.clone();

        info!("The length of the read with barcode is: {}", barcode_read_length);
        info!("The length of the paired read is: {}", paired_read_length);
        (barcode_read, paired_read)
    }
}

#[derive(Getters, Setters, Default, Clone)]
pub struct ReadInfo {
    #[getset(get = "pub", set = "pub")]
    sequence_length: usize,
    #[getset(get = "pub", set = "pub")]
    read_length: usize,
    #[getset(get = "pub", set = "pub")]
    l_position: usize,
    #[getset(get = "pub", set = "pub")]
    flowcell: String,
    #[getset(get = "pub", set = "pub")]
    lane: String,
}

impl ReadInfo {
    pub fn new(
        sequence_length: usize,
        read_length: usize,
        l_position: usize,
        flowcell: String,
        lane: String
    ) -> Self {
        Self { sequence_length, read_length, l_position, flowcell, lane }
    }
}

fn get_lane_from_file(file_name: &PathBuf, loc: usize, sep: char) -> String {
    let tmp: Vec<&str> = file_name.file_name().unwrap().to_str().unwrap().split(sep).collect();
    if tmp.len() > loc {
        let lane = tmp[loc].to_string();
        if lane.starts_with("L0") {
            return lane;
        }
    }
    String::new()
}

fn validate_and_assigne_input_reads(r1: &String, r2: &String) -> (PathBuf, PathBuf, bool) {
    let barcode_read;
    let mut paired_read = String::new();

    if r2.len() == 0 {
        //Single End reads
        barcode_read = r1.to_string();
        info!("Single end read input was detected!");
        info!("Read with Barcode or R1: {}", barcode_read);
    } else {
        paired_read = r1.clone();
        barcode_read = r2.to_string();
        info!("Paired ended read input was detected!");
        info!("Paired read or R1: {}", paired_read);
        info!("Read with Barcode or R2: {}", barcode_read);

        if paired_read.len() == 0 {
            panic!("Input reads are invalid! For single end fastq, use -f or -i parameters!");
        }
        check_file(&paired_read);
    }

    if barcode_read.len() == 0 {
        panic!("Input reads are invalid! check the path {}", barcode_read);
    }
    check_file(&barcode_read);
    (PathBuf::from(paired_read), PathBuf::from(barcode_read), r2.len() == 0)
}

fn get_read_parts(reader: &mut dyn BufRead) -> (String, String, String, String) {
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

fn get_flowcell_lane_info(mgi_header: &String) -> (String, usize, String) {
    let mut l_position = mgi_header.len() - 1;
    for header_chr in mgi_header.chars().rev() {
        if header_chr == 'L' {
            break;
        }

        l_position -= 1;
    }

    if l_position == 0 {
        panic!("Can not find the flowcell id in this header {}!", mgi_header);
    }
    let flowcell = mgi_header[1..l_position].to_string();
    let lane = mgi_header[l_position + 1..l_position + 2].to_string();

    info!("Detected flowcell from the header of the first read is {}.", flowcell);
    (flowcell, l_position, lane)
}

fn get_read_files_from_input_dir(
    in_dir: &String,
    r1_suf: &String,
    r2_suf: &String
) -> (PathBuf, PathBuf, bool) {
    if in_dir.len() > 0 {
        let entries = fs
            ::read_dir(in_dir)
            .expect("can not read directory content!")
            .map(|res| res.map(|e| e.path()))
            .collect::<Result<Vec<_>, io::Error>>()
            .expect("can not collect paths!");
        let mut r1_path: String = String::new();
        let mut r2_path: String = String::new();
        let mut found = 0;
        for path in entries {
            if path.file_name().unwrap().to_str().unwrap().ends_with(r1_suf) {
                r1_path = path.to_str().unwrap().to_string();
                found += 1;
            } else if path.file_name().unwrap().to_str().unwrap().ends_with(r2_suf) {
                r2_path = path.to_str().unwrap().to_string();
                found += 1;
            }
            if found > 1 {
                break;
            }
        }
        if found == 0 {
            panic!(
                "Can not find files that ends with {} or {} under the directory {}",
                in_dir,
                r1_suf,
                r2_suf
            );
        }
        validate_and_assigne_input_reads(&r1_path, &r2_path)
    } else {
        panic!("input directory is not provided!");
    }
}

fn find_info_file<P: AsRef<Path>>(info_file_arg: &P, input_dir: &P, r_file: &P) -> PathBuf {
    if info_file_arg.as_ref().is_file() {
        check_file(info_file_arg);
        return PathBuf::from(info_file_arg.as_ref());
    } else if input_dir.as_ref().is_dir() {
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
    } else if r_file.as_ref().is_file() {
        let tmp_path: PathBuf = r_file.as_ref().with_file_name("BioInfo.csv");
        if tmp_path.exists() {
            return tmp_path;
        }
    }
    return PathBuf::new();
}

fn parse_info_file(info_file_path: &PathBuf) -> (String, String) {
    let mut instrument: String = String::new();
    let mut run: String = String::new();
    if info_file_path.is_file() {
        let file = File::open(info_file_path).unwrap();
        for line in io::BufReader::new(file).lines() {
            if let Ok(inf) = line {
                //info!("{}", inf);
                if inf.starts_with("Machine ID,") {
                    let tmp: Vec<&str> = inf.split(",").collect();
                    instrument = tmp[1].to_string().trim().to_string();
                } else if
                    inf.starts_with("Sequence Date") ||
                    inf.starts_with("Sequence Start Date")
                {
                    let tmp: Vec<&str> = inf.split(",").collect();
                    run = tmp[1].to_string().trim().replace("-", "");
                } else if
                    inf.starts_with("Sequence Time") ||
                    inf.starts_with("Sequence Start Time")
                {
                    let tmp: Vec<&str> = inf.split(",").collect();
                    run.push_str(&tmp[1].to_string().trim().replace(":", ""));
                }
            }
        }
    }
    (instrument, run)
}

fn prepare_output_report_dir(
    ouput_dir_arg: &String,
    report_dir_arg: &String,
    force: bool
) -> (PathBuf, PathBuf) {
    let use_same_dir: bool;
    let output_directory = if ouput_dir_arg.len() == 0 {
        PathBuf::from(&Local::now().format("mgiKit_%Y%m%dT%H%M%S").to_string())
    } else {
        PathBuf::from(ouput_dir_arg)
    };

    let report_directory = if report_dir_arg.len() == 0 {
        info!("The same output directory will be used for reports.");
        use_same_dir = true;
        output_directory.clone()
    } else {
        use_same_dir = false;
        PathBuf::from(report_dir_arg)
    };

    if output_directory.is_dir() {
        if !force {
            panic!(
                "Output directly exists. Use --froce to overwrite their data: {}",
                output_directory.display()
            );
        } else {
            info!(
                "Output directory exists. Data will be overwritten at: {}.",
                output_directory.display()
            );
        }
    } else {
        create_folder(&output_directory);
    }

    if report_directory.is_dir() {
        if !use_same_dir {
            if !force {
                panic!(
                    "Report directly exists. Use --froce to overwrite their data: {}",
                    report_directory.display()
                );
            } else {
                info!(
                    "Report directory exists. Data will be overwritten at: {}.",
                    report_directory.display()
                );
            }
        }
    } else {
        create_folder(&report_directory);
    }

    (output_directory, report_directory)
}
