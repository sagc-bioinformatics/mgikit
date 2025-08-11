use getset::CopyGetters;
use std::{ path::PathBuf, sync::Mutex };
use std::fs::OpenOptions;
use crate::{ RunManager, SampleManager };
use libdeflater::{ Compressor, CompressionLvl };
use crate::file_utils::delete_file;
use log::{ info, warn, debug };
use std::io::Write;

use crate::variables::*;

const MIN_BUFFER_SIZE: usize = 65536;
const MAX_BUFFER_SIZE: usize = 536870912;
const HEADER_TAIL: [u8; 5] = [b':', b'N', b':', b'0', b':'];

pub struct SampleData {
    label: String,
    compressor: Compressor,
    barcode_reads: Option<SampleReads>,
    paired_reads: Option<SampleReads>,
    buffer_info: BufferInfo,
    illumina_header_prefix: String,
}

impl SampleData {
    pub fn new(
        label: String,
        paired_end: bool,
        read2_has_sequence: bool,
        buffer_info: BufferInfo,
        illumina_header_prefix: String
    ) -> Self {
        Self {
            label,
            barcode_reads: match read2_has_sequence {
                false => None,
                true =>
                    Some(
                        SampleReads::new(
                            PathBuf::new(),
                            vec![0; buffer_info.reqiured_output_buffer_size()],
                            vec![0; buffer_info.compression_buffer_size()],
                            0,
                            0
                        )
                    ),
            },
            paired_reads: match paired_end {
                false => None,
                true =>
                    Some(
                        SampleReads::new(
                            PathBuf::new(),
                            vec![0; buffer_info.reqiured_output_buffer_size()],
                            vec![0; buffer_info.compression_buffer_size()],
                            0,
                            0
                        )
                    ),
            },
            compressor: Compressor::new(
                CompressionLvl::new(buffer_info.compression_level() as i32).unwrap()
            ),
            buffer_info,
            illumina_header_prefix,
        }
    }

    pub fn minimal_sample() -> Self {
        Self {
            label: String::new(),
            barcode_reads: None,
            paired_reads: None,
            compressor: Compressor::new(CompressionLvl::new(1_i32).unwrap()),
            buffer_info: BufferInfo::default(),
            illumina_header_prefix: String::new(),
        }
    }

    pub fn is_single_end(&self) -> bool {
        match self.paired_reads {
            None => true,
            Some(_) => false,
        }
    }

    pub fn read2_has_sequence(&self) -> bool {
        match self.barcode_reads {
            None => true,
            Some(_) => false,
        }
    }

    pub fn create_files(
        &mut self,
        lane: &String,
        sample_indx: usize,
        illumina_format: bool,
        output_dir: &PathBuf,
        paired_read_input: bool
    ) {
        let (output_file_r1, output_file_r2) = create_output_file_name(
            &self.label,
            &lane,
            sample_indx,
            illumina_format,
            paired_read_input
        );
        //debug!("create files for sample {}: {}, {} with path {}", sample_indx, lane, illumina_format, output_dir.display());
        match self.paired_reads {
            Some(ref mut sr) => {
                sr.output_file = output_dir.join(output_file_r1);
                //debug!("paired path: {}", sr.output_file.display());
            }
            None => {}
        }
        match self.barcode_reads {
            Some(ref mut sr) => {
                sr.output_file = output_dir.join(output_file_r2);
                //debug!("barcode path: {}", sr.output_file.display());
            }
            None => {}
        };
    }

    pub fn delete_sample_files(&self) {
        match &self.barcode_reads {
            Some(sr) => {
                delete_file(&sr.output_file);
            }
            None => {}
        }
        match &self.paired_reads {
            Some(sr) => {
                delete_file(&sr.output_file);
            }
            None => {}
        };
    }

    pub fn barcode_read_buffer_end(&self) -> usize {
        match &self.barcode_reads {
            Some(sr) => { sr.out_buffer_last }
            None => { 0 }
        }
    }
    pub fn paired_read_buffer_end(&self) -> usize {
        match &self.paired_reads {
            Some(sr) => { sr.out_buffer_last }
            None => { 0 }
        }
    }

    pub fn barcode_read_compression_end(&self) -> usize {
        match &self.barcode_reads {
            Some(sr) => { sr.compression_buffer_last }
            None => { 0 }
        }
    }
    pub fn paired_read_compression_end(&self) -> usize {
        match &self.paired_reads {
            Some(sr) => { sr.compression_buffer_last }
            None => { 0 }
        }
    }

    pub fn print_buffers(&self, start: usize, end: usize, read: usize) {
        if read == 0 || read == 2 {
            match &self.barcode_reads {
                Some(sr) => {
                    sr.print_part(start, end);
                }
                None => {}
            }
        }
        if read == 0 || read == 1 {
            match &self.paired_reads {
                Some(sr) => {
                    sr.print_part(start, end);
                }
                None => {}
            }
        }
    }

    pub fn compress_and_write(&mut self, force: bool, lock: &Mutex<bool>) {
        //debug!("br_comp_end: {}, pr_comp_end: {}, br_out_end: {}, pr_out_end: {}", self.barcode_read_compression_end(), self.paired_read_compression_end(), self.barcode_read_buffer_end(), self.paired_read_buffer_end());
        if
            self.barcode_read_compression_end() >= self.buffer_info.compression_threshold() ||
            self.paired_read_compression_end() >= self.buffer_info.compression_threshold() ||
            force
        {
            match self.barcode_reads {
                Some(ref mut sr) => {
                    sr.compress(&mut self.compressor);
                }
                None => {}
            }
            match self.paired_reads {
                Some(ref mut sr) => {
                    sr.compress(&mut self.compressor);
                }
                None => {}
            };
        }

        if
            self.barcode_read_buffer_end() >= self.buffer_info.writing_threshold() ||
            self.paired_read_buffer_end() >= self.buffer_info.writing_threshold() ||
            force
        {
            let _l = lock.lock().unwrap();
            self.write();
        }
    }

    pub fn write(&mut self) {
        match self.barcode_reads {
            Some(ref mut sr) => {
                sr.write();
            }
            None => {}
        }
        match self.paired_reads {
            Some(ref mut sr) => {
                sr.write();
            }
            None => {}
        };
    }

    pub fn add_barcode_reads(&mut self, data: &[u8]) {
        match self.barcode_reads {
            Some(ref mut sr) => {
                sr.add_reads(data);
            }
            None => {}
        };
    }

    pub fn add_paired_reads(&mut self, data: &[u8]) {
        match self.paired_reads {
            Some(ref mut sr) => {
                sr.add_reads(data);
            }
            None => {}
        };
    }

    pub fn write_illumina_header(
        &mut self,
        mgi_read_header: &[u8],
        umi: &[u8],
        l_position: usize,
        sep_position: usize,
        sb_header: bool,
        barcode_read: bool,
        paired_read: bool
    ) {
        if barcode_read {
            match self.barcode_reads {
                Some(ref mut sr) => {
                    //debug!("Barcode reads - compre buffer: exp size {} - last {} - actual size: {}", self.buffer_info.compression_buffer_size, sr.compression_buffer_last, sr.compression_buffer.len());
                    sr.add_reads(&self.illumina_header_prefix.as_bytes());
                    sr.compression_buffer_last = write_illumina_header(
                        &mut sr.compression_buffer,
                        sr.compression_buffer_last,
                        mgi_read_header,
                        umi,
                        l_position,
                        sep_position,
                        sb_header
                    );
                }
                None => {}
            };
        }
        if paired_read {
            match self.paired_reads {
                Some(ref mut sr) => {
                    //debug!("Paired reads - compre buffer: exp size {} - last {} - actual size: {}", self.buffer_info.compression_buffer_size, sr.compression_buffer_last, sr.compression_buffer.len());
                    sr.add_reads(&self.illumina_header_prefix.as_bytes());
                    sr.compression_buffer_last = write_illumina_header(
                        &mut sr.compression_buffer,
                        sr.compression_buffer_last,
                        mgi_read_header,
                        umi,
                        l_position,
                        sep_position,
                        sb_header
                    );
                }
                None => {}
            };
        }
    }

    pub fn get_barcode_read_compression_end(&self) -> usize {
        match &self.barcode_reads {
            Some(sr) => { sr.compression_buffer_last }
            None => { usize::MAX }
        }
    }
    pub fn get_paired_read_compression_end(&self) -> usize {
        match &self.paired_reads {
            Some(sr) => { sr.compression_buffer_last }
            None => { usize::MAX }
        }
    }
    pub fn fix_paired_read_header(&mut self, tail_offset: usize, val: u8) {
        match self.paired_reads {
            Some(ref mut sr) => {
                sr.compression_buffer[sr.compression_buffer_last - tail_offset] = val;
            }
            None => {}
        };
    }

    pub fn copy_header_2_paired(&mut self, start: usize) {
        match self.paired_reads {
            Some(ref mut pr) => {
                match self.barcode_reads {
                    Some(ref mut br) => {
                        pr.compression_buffer[
                            pr.compression_buffer_last..pr.compression_buffer_last +
                                br.compression_buffer_last -
                                start
                        ].copy_from_slice(
                            &br.compression_buffer[start..br.compression_buffer_last]
                        );
                        pr.compression_buffer_last += br.compression_buffer_last - start;
                    }
                    None => {}
                };
            }
            None => {}
        };
    }
}

pub struct SampleReads {
    pub output_file: PathBuf,
    pub out_buffer: Vec<u8>,
    pub compression_buffer: Vec<u8>,
    pub out_buffer_last: usize,
    pub compression_buffer_last: usize,
}

impl SampleReads {
    pub fn new(
        output_file: PathBuf,
        out_buffer: Vec<u8>,
        compression_buffer: Vec<u8>,
        out_buffer_last: usize,
        compression_buffer_last: usize
    ) -> Self {
        Self {
            output_file,
            out_buffer,
            compression_buffer,
            out_buffer_last,
            compression_buffer_last,
        }
    }

    /* 
    pub fn compress_and_write(&mut self, compressor: &mut Compressor, compression_threshold: usize, writing_threshold: usize, force: bool){
        if self.compression_buffer_last >= compression_threshold || force {
            self.out_buffer_last = compress_buffer(&self.compression_buffer, self.compression_buffer_last, &mut self.out_buffer, self.out_buffer_last, compressor);
            self.compression_buffer_last = 0;
        }
        if self.out_buffer_last >= writing_threshold || force {
            write_data(&self.out_buffer, self.out_buffer_last, &self.output_file);
            self.out_buffer_last = 0;
        }       
    }
    */

    pub fn print_part(&self, start: usize, end: usize) {
        println!(
            "compression_buffer: *{}*",
            String::from_utf8(self.compression_buffer[start..end].to_vec()).unwrap()
        )
    }

    pub fn compress(&mut self, compressor: &mut Compressor) {
        self.out_buffer_last = compress_buffer(
            &self.compression_buffer,
            self.compression_buffer_last,
            &mut self.out_buffer,
            self.out_buffer_last,
            compressor
        );
        self.compression_buffer_last = 0;
    }

    pub fn write(&mut self) {
        write_data(&self.out_buffer, self.out_buffer_last, &self.output_file);
        self.out_buffer_last = 0;
    }

    pub fn add_reads(&mut self, data: &[u8]) {
        //debug!("adding read with length {}", data.len());
        self.compression_buffer[
            self.compression_buffer_last..self.compression_buffer_last + data.len()
        ].copy_from_slice(data);
        self.compression_buffer_last += data.len();
    }
}

#[derive(CopyGetters, Clone, Default)]
#[getset(get_copy = "pub")]
pub struct BufferInfo {
    writing_buffer_size: usize,
    compression_buffer_size: usize,
    compression_level: u32,
    reqiured_output_buffer_size: usize,
    compression_threshold: usize,
    writing_threshold: usize,
}

impl BufferInfo {
    pub fn new(
        writing_buffer_size: usize,
        compression_buffer_size: usize,
        compression_level: u32,
        compression_threshold: usize,
        writing_threshold: usize
    ) -> Self {
        info!("Output buffer size: {}", writing_buffer_size);
        info!("Compression buffer size: {}", compression_buffer_size);
        info!("Compression level: {}. (0 no compression but fast, 12 best compression but slow.)", compression_level);
        let mut compressor = Compressor::new(
            CompressionLvl::new(compression_level as i32).unwrap()
        );
        if compression_buffer_size > writing_buffer_size {
            panic!(
                "Compression buffer size '--compression-buffer-size' should be less than Writing buffer size ('--writing-buffer-size')."
            );
        }
        if writing_buffer_size < MIN_BUFFER_SIZE {
            panic!("Writing buffer size '--writing-buffer-size' will be increased to minimal allowed value ({}).", MIN_BUFFER_SIZE);
        } else if writing_buffer_size > MAX_BUFFER_SIZE {
            panic!("Writing buffer size '--writing-buffer-size' will be reduced to the maximum allowed value ({}).", MAX_BUFFER_SIZE);
        }

        let reqiured_output_buffer_size =
            writing_buffer_size + compressor.gzip_compress_bound(compression_buffer_size);
        debug!("Compression buffer flush threshold: {}", compression_threshold);
        debug!("Output buffer flush threshold: {}", writing_threshold);

        Self {
            writing_buffer_size: writing_buffer_size,
            compression_buffer_size,
            compression_level,
            reqiured_output_buffer_size,
            compression_threshold,
            writing_threshold,
        }
    }

    pub fn calculate_final_writing_buffer_size(&mut self, max_buffer_size: usize) {
        if self.writing_buffer_size > max_buffer_size {
            warn!("Writing buffer size has been reduced to {} as there is no enough memory!", max_buffer_size);
            self.reqiured_output_buffer_size =
                max_buffer_size +
                Compressor::new(
                    CompressionLvl::new(self.compression_level as i32).unwrap()
                ).gzip_compress_bound(self.compression_buffer_size);

            self.writing_threshold = max_buffer_size;
        }
    }
}

fn compress_buffer(
    source: &Vec<u8>,
    source_size: usize,
    destination: &mut Vec<u8>,
    destination_size: usize,
    compressor: &mut Compressor
) -> usize {
    //debug!("Compressing {} bytes, on {}", source_size, destination_size);
    if source_size > 0 {
        return
            destination_size +
            compressor
                .gzip_compress(&source[0..source_size], &mut destination[destination_size..])
                .unwrap();
    }
    return destination_size;
}

fn write_data(out_buffer: &Vec<u8>, buffer_size: usize, output_file_path: &PathBuf) {
    //debug!("writing {} bytes into {}", buffer_size, output_file_path.display());
    if buffer_size > 0 {
        let mut curr_writer = OpenOptions::new()
            .append(true)
            .create(true)
            .open(output_file_path)
            .expect("couldn't create output");
        curr_writer.write_all(&out_buffer[..buffer_size]).unwrap();
        curr_writer.flush().expect("couldn't flush output");
    }
}

pub fn calculate_reqiured_memory(
    total_samples: usize,
    writing_buffer_size: usize,
    compression_buffer_size: usize,
    single_read_input: bool
) -> f64 {
    if single_read_input {
        (total_samples * (writing_buffer_size + 2 * compression_buffer_size)) as f64
    } else {
        2_f64 * ((total_samples * (writing_buffer_size + 2 * compression_buffer_size)) as f64)
    }
}

pub fn calculate_largest_buffer_size(
    available_memory: f64,
    mut total_samples: usize,
    compression_buffer_size: usize,
    single_end: bool,
    processors: usize
) -> usize {
    if single_end {
        total_samples *= 2;
    }
    (2_usize).pow(
        (
            available_memory / ((total_samples * processors) as f64) -
            2_f64 * (compression_buffer_size as f64)
        )
            .log2()
            .floor() as u32
    )
}

fn create_output_file_name(
    sample_name: &String,
    lane: &String,
    sample_index: usize,
    illumina_format: bool,
    paired_read_input: bool
) -> (String, String) {

    let br_suff = if paired_read_input{"R2"}else{"R1"};

    if illumina_format {
        if sample_index == usize::MAX {
            return (
                format!("{}_{}_R1_001.fastq.gz", sample_name, lane),
                format!("{}_{}_{}_001.fastq.gz", sample_name, lane, br_suff),
            );
        } else {
            return (
                format!("{}_S{}_{}_R1_001.fastq.gz", sample_name, sample_index, lane),
                format!("{}_S{}_{}_{}_001.fastq.gz", sample_name, sample_index, lane, br_suff),
            );
        }
    } else {
        return (
            format!("{}_{}_R1.fastq.gz", sample_name, lane),
            format!("{}_{}_{}.fastq.gz", sample_name, lane, br_suff),
        );
        //return (format!("{}_R1.fastq.gz", sample_name), format!("{}_R2.fastq.gz", sample_name));
    }
}

fn write_illumina_header(
    output_buffer: &mut [u8],
    mut buffer_end: usize,
    mgi_read_header: &[u8],
    umi: &[u8],
    l_position: usize,
    sep_position: usize,
    sb_header: bool
) -> usize {
    output_buffer[buffer_end] = mgi_read_header[l_position + 1];
    buffer_end += 1;

    output_buffer[buffer_end] = b':';
    buffer_end += 1;

    for i in l_position + 10..sep_position {
        if mgi_read_header[i] != b'0' {
            output_buffer[buffer_end..buffer_end + sep_position - i].copy_from_slice(
                &mgi_read_header[i..sep_position]
            );
            buffer_end += sep_position - i;
            break;
        }
    }

    output_buffer[buffer_end] = b':';
    buffer_end += 1;

    if mgi_read_header[l_position + 3] != b'0' {
        output_buffer[buffer_end] = mgi_read_header[l_position + 3];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 4];
        output_buffer[buffer_end + 2] = mgi_read_header[l_position + 5];
        buffer_end += 3;
    } else if mgi_read_header[l_position + 4] != b'0' {
        output_buffer[buffer_end] = mgi_read_header[l_position + 4];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 5];
        buffer_end += 2;
    } else {
        output_buffer[buffer_end] = mgi_read_header[l_position + 5];
        buffer_end += 1;
    }

    output_buffer[buffer_end] = b':';
    buffer_end += 1;

    //debug!("")
    if mgi_read_header[l_position + 7] != b'0' {
        output_buffer[buffer_end] = mgi_read_header[l_position + 7];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 8];
        output_buffer[buffer_end + 2] = mgi_read_header[l_position + 9];
        buffer_end += 3;
    } else if mgi_read_header[l_position + 8] != b'0' {
        output_buffer[buffer_end] = mgi_read_header[l_position + 8];
        output_buffer[buffer_end + 1] = mgi_read_header[l_position + 9];
        buffer_end += 2;
    } else {
        output_buffer[buffer_end] = mgi_read_header[l_position + 9];
        buffer_end += 1;
    }

    if sb_header {
        if umi.len() > 0 {
            output_buffer[buffer_end] = b':';
            buffer_end += 1;
            output_buffer[buffer_end..buffer_end + umi.len()].copy_from_slice(&umi[..]);
            buffer_end += umi.len();
        }
        output_buffer[
            buffer_end..buffer_end + mgi_read_header.len() - sep_position - 1
        ].copy_from_slice(&mgi_read_header[sep_position..mgi_read_header.len() - 1]);

        buffer_end += mgi_read_header.len() - sep_position - 1;
    } else {
        if umi.len() > 0 {
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

pub fn create_sample_data_list(
    sample_manager: &SampleManager,
    run_manager: &RunManager,
    buffer_info: &BufferInfo,
    read2_has_sequence: bool,
    illumina_format: bool
) -> Vec<SampleData> {
    let mut sample_data_list: Vec<SampleData> = Vec::new();
    let writing_samples = sample_manager.writing_samples();
    let unique_samples_ids = sample_manager.unique_samples_ids();
    let sample_information: &Vec<Vec<String>> = sample_manager.sample_information();
    let total_samples = sample_manager.get_sample_count();
    let undetermined_label_id = total_samples - 2;
    //debug!("creating samples: {} - {} - {} - {}", total_samples, read2_has_sequence, illumina_format, run_manager.output_dir().display());
    let illumina_header = if illumina_format {
        run_manager.create_illumina_header_prefix()
    } else {
        String::new()
    };
    {
        let mut sample_data: SampleData;
        for i in 0..total_samples {
            if writing_samples[i] == i {
                //debug!("creating complete samples: {} - {}", writing_samples[i], i);
                sample_data = SampleData::new(
                    sample_information[i][SAMPLE_COLUMN].clone(),
                    run_manager.paired_read_input(),
                    read2_has_sequence || i >= undetermined_label_id,
                    buffer_info.clone(),
                    illumina_header.clone()
                );
                if i < undetermined_label_id && illumina_format {
                    sample_data.create_files(
                        run_manager.lane(),
                        unique_samples_ids[i] + 1,
                        true,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                } else if i >= undetermined_label_id && illumina_format {
                    sample_data.create_files(
                        run_manager.lane(),
                        usize::MAX,
                        true,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                } else {
                    sample_data.create_files(
                        run_manager.lane(),
                        unique_samples_ids[i] + 1,
                        false,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                }
            } else {
                //debug!("creating empty samples: {} - {}", writing_samples[i], i);
                sample_data = SampleData::new(
                    sample_information[i][SAMPLE_COLUMN].clone(),
                    false,
                    false,
                    buffer_info.clone(),
                    illumina_header.clone()
                );
            }
            sample_data_list.push(sample_data);
        }
    }
    sample_data_list
}

pub fn clean_output_directory(
    sample_manager: &SampleManager,
    run_manager: &RunManager,
    buffer_info: &BufferInfo,
    read2_has_sequence: bool,
    illumina_format: bool
) {
    let writing_samples = sample_manager.writing_samples();
    let unique_samples_ids = sample_manager.unique_samples_ids();
    let sample_information: &Vec<Vec<String>> = sample_manager.sample_information();
    let total_samples = sample_manager.get_sample_count();
    let undetermined_label_id = total_samples - 2;
    let illumina_header = if illumina_format {
        run_manager.create_illumina_header_prefix()
    } else {
        String::new()
    };
    {
        let mut sample_data: SampleData;
        for i in 0..total_samples {
            if writing_samples[i] == i {
                sample_data = SampleData::new(
                    sample_information[i][SAMPLE_COLUMN].clone(),
                    run_manager.paired_read_input(),
                    read2_has_sequence || i >= undetermined_label_id,
                    buffer_info.clone(),
                    illumina_header.clone()
                );
                if i < undetermined_label_id && illumina_format {
                    sample_data.create_files(
                        run_manager.lane(),
                        unique_samples_ids[i] + 1,
                        true,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                    sample_data.delete_sample_files();
                } else if i >= undetermined_label_id && illumina_format {
                    sample_data.create_files(
                        run_manager.lane(),
                        usize::MAX,
                        true,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                    sample_data.delete_sample_files();
                } else {
                    sample_data.create_files(
                        run_manager.lane(),
                        unique_samples_ids[i] + 1,
                        false,
                        run_manager.output_dir(),
                        run_manager.paired_read_input()
                    );
                    sample_data.delete_sample_files();
                }
            }
        }
    }
}
