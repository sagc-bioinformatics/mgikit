use std::fs;
use std::path::{ Path, PathBuf };
use log::info;
use core::panic;
use flate2::read::MultiGzDecoder;
use flate2::{ Decompress, FlushDecompress, Status };
use std::io::{ self, BufReader, BufWriter, Read, Write };
use std::fs::File;
use std::thread::{ self, JoinHandle };
use crossbeam_channel::{ Receiver, Sender };
use memchr::memchr_iter;
use log::debug;
use std::mem;

pub struct RawReader {
    sender: Sender<(usize, Vec<u8>)>,
    receiver: Receiver<(usize, Vec<u8>)>,
    buffer: Vec<u8>,
    offset: usize,
    end: usize,
    done: bool,
}

impl RawReader {
    pub fn new(
        receiver: Receiver<(usize, Vec<u8>)>,
        sender: Sender<(usize, Vec<u8>)>,
        buffer_size: usize
    ) -> Self {
        Self {
            sender,
            receiver,
            buffer: vec![0_u8; buffer_size],
            offset: 0,
            end: 0,
            done: false,
        }
    }
}

impl Read for RawReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.offset == self.end && !self.done {
            match self.receiver.recv() {
                Ok((bytes, mut chunk)) => {
                    mem::swap(&mut self.buffer, &mut chunk);
                    match self.sender.send((0, chunk)) {
                        Ok(_) => {
                            //debug!("Sending: {} bytes", last_sent_endline);
                            //debug!("kept: {}",extra_len);
                        }
                        Err(e) => println!("error: {e:?}"),
                    }
                    self.offset = 0;
                    self.end = bytes;
                }
                Err(_) => {
                    panic!("Error receiving!");
                }
            }
        }

        if self.offset < self.end {
            let n = (self.end - self.offset).min(out.len());
            out[..n].copy_from_slice(&self.buffer[self.offset..self.offset + n]);
            self.offset += n;
            Ok(n)
        } else {
            Ok(0) // EOF
        }
    }
}

fn read_bytes_in_reads<R: Read>(
    reader: &mut R,
    buffer: &mut [u8],
    minimum: usize,
    last_byte: &mut usize
) -> (usize, Vec<usize>) {
    let mut curr_bytes: usize = 0;
    let mut lines: Vec<usize> = Vec::with_capacity(minimum * 2);
    lines.extend(memchr_iter(b'\n', &buffer[..*last_byte]));
        
    loop {
        debug!(
            "buffer length: {}, starting from {}, current read {}",
            buffer.len(),
            last_byte,
            curr_bytes
        );
        curr_bytes = reader.read(&mut buffer[*last_byte..]).unwrap();
        lines.extend(memchr_iter(b'\n', &buffer[*last_byte..*last_byte + curr_bytes]).map(|it| it + *last_byte));
        *last_byte += curr_bytes;
        //println!(" reads: {} - lines: {}  - curr: {}", minimum * 4, lines.len(), curr_bytes);
        if curr_bytes == 0 || lines.len() >= minimum * 4 {
            //debug!("total lines: {} - no more", line_cnt);
            break;
        }
    }
    //let last_10 = &lines[lines.len().saturating_sub(10)..];
    //println!("{}  - {:?}", lines.len(), last_10);
    //lines = Vec::new();
    //lines.extend(memchr_iter(b'\n', &buffer[..*last_byte]).collect::<Vec<usize>>());
    //let last_10 = &lines[lines.len().saturating_sub(10)..];
    //println!("{}  - {:?}", lines.len(), last_10);
    
    //debug!("total lines: {} - still more", line_cnt);
    return (*last_byte, lines);
}

pub fn fill_send_buffers<R: Read>(
    full_sender: &Sender<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver: &Receiver<(usize, Vec<u8>)>,
    empty_sender: &Sender<(usize, Vec<u8>)>,
    reader: &mut R,
    read_cnt: usize,
    extra: &mut Vec<u8>,
    extra_len: &mut usize,
    force: bool
) -> (usize, usize) {
    let mut total_bytes = 0;
    match empty_receiver.recv() {
        Ok((_, mut buffer)) => {
            if *extra_len > 0 {
                buffer[..*extra_len].copy_from_slice(&extra[..*extra_len]);
                total_bytes = *extra_len;
                *extra_len = 0;
            }

            let (last_read_bytes, mut lines) = read_bytes_in_reads(
                reader,
                &mut buffer[..],
                read_cnt,
                &mut total_bytes
            );
            //let mut lines = memchr_iter(b'\n', &buffer[..total_bytes]);
            //debug!("total read bytes: {}  -  lines: {}  - needed lines: {}", total_bytes, lines.count(), read_cnt * 4);
            //let mut lines = memchr_iter(b'\n', &buffer[..total_bytes]);
            if lines.len() >= read_cnt * 4 || force {
                if lines.len() >= read_cnt * 4 {
                    lines.truncate(read_cnt * 4);
                    *extra_len = total_bytes - lines.last().unwrap() - 1;
                    extra[..*extra_len].copy_from_slice(
                        &buffer[lines.last().unwrap() + 1..total_bytes]
                    );
                    total_bytes = lines.last().unwrap() + 1;
                }
                //println!("reader lines: {}", lines.count());
                match full_sender.send((total_bytes, buffer, lines)) {
                    Ok(_) => {
                        //debug!("Sending: {} bytes", last_sent_endline);
                        //debug!("kept: {}",extra_len);
                    }
                    Err(e) => println!("error: {e:?}"),
                }
                return (total_bytes, last_read_bytes);
            } else {
                *extra_len = total_bytes;
                extra[..*extra_len].copy_from_slice(&buffer[..total_bytes]);
                match empty_sender.send((0, buffer)) {
                    Ok(_) => {
                        //debug!("Sending: {} bytes", last_sent_endline);
                        //debug!("kept: {}",extra_len);

                    }
                    Err(e) => println!("error: {e:?}"),
                }
                return (0, last_read_bytes);
            }
        }
        Err(_) => {
            debug!("Error while sending!");
            panic!("error!")
        }
    }
}

pub fn send_raw_data_buffers<R: Read>(
    full_sender: &Sender<(usize, Vec<u8>)>,
    empty_receiver: &Receiver<(usize, Vec<u8>)>,
    reader: &mut R
) -> bool {
    let mut total_bytes = 0;
    match empty_receiver.recv() {
        Ok((_, mut buffer)) => {
            let _ = read_bytes(reader, &mut buffer[..], &mut total_bytes);
            match full_sender.send((total_bytes, buffer)) {
                Ok(_) => {
                    //debug!("Sending: {} bytes", last_sent_endline);
                    //debug!("kept: {}",extra_len);

                }
                Err(e) => println!("error: {e:?}"),
            }
            if total_bytes == 0 {
                return false;
            }
        }
        Err(_) => {
            debug!("Error while sending!");
        }
    }
    return true;
}

pub fn check_file<P: AsRef<Path>>(path: &P) {
    if !path.as_ref().is_file() {
        panic!("File is not accessible: {}", path.as_ref().display());
    }
}

pub fn create_folder<P: AsRef<Path>>(path: &P) {
    info!("A directory is created: {}", path.as_ref().display());
    fs::create_dir_all(path.as_ref()).unwrap();
}

pub fn get_buf_reader(input_file: &PathBuf) -> BufReader<MultiGzDecoder<File>> {
    BufReader::new(MultiGzDecoder::new(File::open(input_file).expect("Could not open the file")))
}

/*
pub fn get_reader(input_file: &PathBuf) -> Box<dyn Read> {
    let (reader, _) = niffler
        ::get_reader(Box::new(std::fs::File::open(input_file).unwrap()))
        .unwrap();
    reader
}
*/

pub fn delete_file<P: AsRef<Path>>(path: &P) {
    if path.as_ref().exists() {
        fs::remove_file(path).unwrap();
    }
}

pub fn create_output_file<P: AsRef<Path>>(file_path: &P) -> File {
    File::create(file_path.as_ref()).expect("couldn't create output")
}
pub fn get_buf_writer<P: AsRef<Path>>(file_path: &P) -> BufWriter<File> {
    BufWriter::new(File::create(file_path.as_ref()).expect("couldn't create output"))
}

pub fn get_gzip_reader<P: AsRef<Path>>(file_path: &P) -> MultiGzDecoder<File> {
    flate2::read::MultiGzDecoder::new(std::fs::File::open(file_path).unwrap())
}

pub fn get_reader<P: AsRef<Path>>(file_path: &P) -> File {
    std::fs::File::open(file_path).unwrap()
}

pub fn parallel_reader_thread(
    paired_reads: PathBuf,
    barcode_reads: PathBuf,
    batch_size: usize,
    read_rb: bool,
    read_rp: bool,
    full_sender_rb: Sender<(usize, Vec<u8>, Vec<usize>)>,
    full_sender_rp: Sender<(usize, Vec<u8>, Vec<usize>)>,
    empty_sender_rb: Sender<(usize, Vec<u8>)>,
    empty_sender_rp: Sender<(usize, Vec<u8>)>,
    full_receiver_rb: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_receiver_rp: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver_rb: Receiver<(usize, Vec<u8>)>,
    empty_receiver_rp: Receiver<(usize, Vec<u8>)>,
    full_sender: Sender<(usize, Vec<u8>, Vec<usize>, usize, Vec<u8>, Vec<usize>)>,
    processing_threads: usize,
    buffer_size: usize,
    paired_input: bool,
    main_sender: bool
) -> JoinHandle<()> {
    thread::spawn(move || {
        info!("Reader thread has started!");
        let mut reader_barcode_read = if read_rb {
            Some(get_gzip_reader(&barcode_reads))
        } else {
            None
        };
        let mut reader_paired_read = if read_rp {
            Some(get_gzip_reader(&paired_reads))
        } else {
            None
        };

        let mut extra_rb: Vec<u8> = if read_rb { vec![b'0'; buffer_size] } else { Vec::new() };
        let mut extra_len_rb: usize = 0;

        let mut extra_rp: Vec<u8> = if read_rp { vec![b'0'; buffer_size] } else { Vec::new() };
        let mut extra_len_rp: usize = 0;

        let mut readers_finished = false;
        loop {
            //debug!("--------------------------------");
            if read_rb {
                match reader_barcode_read {
                    Some(ref mut reader_barcode) => {
                        let (sent_bytes, _) = fill_send_buffers(
                            &full_sender_rb,
                            &empty_receiver_rb,
                            &empty_sender_rb,
                            reader_barcode,
                            batch_size,
                            &mut extra_rb,
                            &mut extra_len_rb,
                            true
                        );
                        if !main_sender && sent_bytes == 0 {
                            readers_finished = true;
                        }
                        debug!("Sent rb full buffer!");
                    }
                    None => {}
                }
            }

            if read_rp {
                match reader_paired_read {
                    Some(ref mut reader_paired) => {
                        let (sent_bytes, _) = fill_send_buffers(
                            &full_sender_rp,
                            &empty_receiver_rp,
                            &empty_sender_rp,
                            reader_paired,
                            batch_size,
                            &mut extra_rp,
                            &mut extra_len_rp,
                            true
                        );
                        if !main_sender && sent_bytes == 0 {
                            readers_finished = true;
                        }
                        debug!("Sent rp full buffer!");
                    }
                    None => {}
                }
            }

            if main_sender {
                let (read_bytes2, buffer2, lines_rb) = full_receiver_rb.recv().unwrap();
                let (read_bytes1, buffer1, lines_rp) = if paired_input {
                    full_receiver_rp.recv().unwrap()
                } else {
                    (0, Vec::new(), Vec::new())
                };
                match
                    full_sender.send((
                        read_bytes2,
                        buffer2,
                        lines_rb,
                        read_bytes1,
                        buffer1,
                        lines_rp,
                    ))
                {
                    Ok(_) => {
                        debug!("Sending full {} - {}", read_bytes1, read_bytes2);
                    } //println!("All good: {v:?}"),
                    Err(e) => println!("error: {e:?}"),
                }
                if (paired_input && read_bytes1 == 0) || read_bytes2 == 0 {
                    readers_finished = true;
                }
            }
            //println!("reading: {}   ---  {}", read_rb, read_rp);
            //panic!(" ------------");
            if readers_finished {
                if main_sender {
                    for i in 1..processing_threads {
                        match
                            full_sender.send((0, Vec::new(), Vec::new(), 0, Vec::new(), Vec::new()))
                        {
                            Ok(_) => {
                                debug!("Sending finish signal {}", i);
                            } //println!("All good: {v:?}"),
                            Err(e) => println!("error: {e:?}"),
                        };
                    }
                    info!("Readers threads are done!");
                } else {
                    info!("Secondary reader thread is done!");
                }
                break;
            }
        }
    })
}

pub fn read_bytes<R: Read>(reader: &mut R, buffer: &mut [u8], last_byte: &mut usize) -> bool {
    let mut curr_bytes: usize;
    loop {
        curr_bytes = reader.read(&mut buffer[*last_byte..]).unwrap();
        *last_byte += curr_bytes;

        if *last_byte == buffer.len() || curr_bytes == 0 {
            if curr_bytes == 0 {
                return false;
            }
            return true;
        }
    }
}

pub fn read_buffers<R: Read>(
    mut bytes: usize,
    buffer: &mut Vec<u8>,
    reader_op: &mut Option<R>
) -> usize {
    if bytes >= 10000 {
        return bytes;
    }
    match reader_op {
        Some(ref mut reader) => {
            read_bytes(reader, buffer, &mut bytes);
        }
        None => {
            panic!("There should be a reader here!");
        }
    }

    //let lines_itr = memchr_iter(b'\n', &buffer[..bytes]);
    bytes
}

pub fn write_file<P: AsRef<Path>>(file_path: &P, content: &String) {
    File::create(file_path.as_ref())
        .expect("couldn't create output")
        .write_all(content.as_bytes())
        .unwrap();
}
/*
pub fn parallel_reader_decompressor_thread(
    reads_path: PathBuf,
    batch_size: usize,
    empty_raw_sender: Sender<(usize, Vec<u8>)>,
    empty_raw_receiver: Receiver<(usize, Vec<u8>)>,
    full_raw_sender: Sender<(usize, Vec<u8>)>,
    full_raw_receiver: Receiver<(usize, Vec<u8>)>,
    full_receiver_rb: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_receiver_rp: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_sender: Sender<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver: Receiver<(usize, Vec<u8>)>,
    empty_sender: Sender<(usize, Vec<u8>)>,
    full_sender_paired: Sender<(usize, Vec<u8>, Vec<usize>, usize, Vec<u8>, Vec<usize>)>,
    processing_threads: usize,
    buffer_size: usize,
    paired_input: bool,
    main_sender: bool
) -> JoinHandle<()> {
    thread::spawn(move || {
        info!("Reader thread has started!");
        let decoder_thread = thread::spawn(move || {
            let mut decoder = Decompress::new(true);
            //let mut read_bytes ;
            //let mut buffer = Vec::new();
            let mut extra_len = 0;
            let mut extra = vec![0_u8; buffer_size];
            let mut raw_start;
            let mut decompressed_bytes;
            let mut sent_buffer;
            let mut reset_decompressor;
            loop{
                let (read_bytes, buffer) = full_raw_receiver.recv().unwrap();
                raw_start = 0;
                if read_bytes == 0{
                    send_buffer(
                        &full_sender,
                        &empty_receiver,
                        &extra[..extra_len]
                    );
                    break;
                }
                loop{
                    (decompressed_bytes, sent_buffer, reset_decompressor) = decode_send_buffer(
                        &full_sender,
                        &empty_receiver,
                        &empty_sender,
                        &mut decoder,
                        &buffer[raw_start..read_bytes],
                        batch_size,
                        &mut extra,
                        &mut extra_len,
                        read_bytes == 0
                    );
                    debug!("decompressed_bytes: {} - sent_buffer: {} - reset_decompressor: {}", decompressed_bytes, sent_buffer, reset_decompressor);
                    if reset_decompressor{
                        decoder = Decompress::new(true);
                    }
                    raw_start += decompressed_bytes;
                    
                    if main_sender && sent_buffer {
                        let (read_bytes2, buffer2, lines_rb) = full_receiver_rb.recv().unwrap();
                        let (read_bytes1, buffer1, lines_rp) = if paired_input {
                            full_receiver_rp.recv().unwrap()
                        } else {
                            (0, Vec::new(), Vec::new())
                        };
                        match
                            full_sender_paired.send((
                                read_bytes2,
                                buffer2,
                                lines_rb,
                                read_bytes1,
                                buffer1,
                                lines_rp,
                            ))
                        {
                            Ok(_) => {
                                debug!("Sending full {} - {}", read_bytes1, read_bytes2);
                            } //println!("All good: {v:?}"),
                            Err(e) => println!("error: {e:?}"),
                        }
                    }
                    if raw_start == read_bytes
                    {
                        break;
                    }
                }
                match empty_raw_sender.send((0, buffer)) {
                    Ok(_) => {
                        //debug!("Sending: {} bytes", last_sent_endline);
                        //debug!("kept: {}",extra_len);

                    }
                    Err(e) => println!("error: {e:?}"),
                }
            }
            if main_sender {
                for i in 1..processing_threads {
                    match
                        full_sender_paired.send((0, Vec::new(), Vec::new(), 0, Vec::new(), Vec::new()))
                    {
                        Ok(_) => {
                            debug!("Sending finish signal {}", i);
                        } //println!("All good: {v:?}"),
                        Err(e) => println!("error: {e:?}"),
                    };
                }
                info!("Readers threads are done!");
            } else {
                info!("Secondary reader thread is done!");
            }
        });
        let mut data_reader = get_reader(&reads_path);
        loop {
            if ! send_raw_data_buffers(
                &full_raw_sender,
                &empty_raw_receiver,
                &mut data_reader,
            ){break;}
        }
        decoder_thread.join().unwrap();
    })
}
*/
pub fn parallel_reader_decompressor_thread(
    reads_path: PathBuf,
    batch_size: usize,
    empty_raw_sender: Sender<(usize, Vec<u8>)>,
    empty_raw_receiver: Receiver<(usize, Vec<u8>)>,
    full_raw_sender: Sender<(usize, Vec<u8>)>,
    full_raw_receiver: Receiver<(usize, Vec<u8>)>,
    full_receiver_rb: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_receiver_rp: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_sender: Sender<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver: Receiver<(usize, Vec<u8>)>,
    empty_sender: Sender<(usize, Vec<u8>)>,
    full_sender_paired: Sender<(usize, Vec<u8>, Vec<usize>, usize, Vec<u8>, Vec<usize>)>,
    processing_threads: usize,
    buffer_size: usize,
    paired_input: bool,
    main_sender: bool
) -> JoinHandle<()> {
    thread::spawn(move || {
        info!("Reader thread has started!");
        let decoder_thread = thread::spawn(move || {
            let mut extra_len = 0;
            let mut extra = vec![0_u8; buffer_size];
            let mut decoder = MultiGzDecoder::new(
                RawReader::new(full_raw_receiver, empty_raw_sender, buffer_size)
            );
            loop {
                let (sent_bytes, _) = fill_send_buffers(
                    &full_sender,
                    &empty_receiver,
                    &empty_sender,
                    &mut decoder,
                    batch_size,
                    &mut extra,
                    &mut extra_len,
                    true
                );

                if main_sender {
                    let (read_bytes2, buffer2, lines_rb) = full_receiver_rb.recv().unwrap();
                    let (read_bytes1, buffer1, lines_rp) = if paired_input {
                        full_receiver_rp.recv().unwrap()
                    } else {
                        (0, Vec::new(), Vec::new())
                    };
                    match
                        full_sender_paired.send((
                            read_bytes2,
                            buffer2,
                            lines_rb,
                            read_bytes1,
                            buffer1,
                            lines_rp,
                        ))
                    {
                        Ok(_) => {
                            debug!("Sending full {} - {}", read_bytes1, read_bytes2);
                        } //println!("All good: {v:?}"),
                        Err(e) => println!("error: {e:?}"),
                    }
                }
                if sent_bytes == 0 {
                    break;
                }
            }
            if main_sender {
                for i in 1..processing_threads {
                    match
                        full_sender_paired.send((
                            0,
                            Vec::new(),
                            Vec::new(),
                            0,
                            Vec::new(),
                            Vec::new(),
                        ))
                    {
                        Ok(_) => {
                            debug!("Sending finish signal {}", i);
                        } //println!("All good: {v:?}"),
                        Err(e) => println!("error: {e:?}"),
                    };
                }
                info!("Readers threads are done!");
            } else {
                info!("Secondary reader thread is done!");
            }
        });
        let mut data_reader = get_reader(&reads_path);
        loop {
            if !send_raw_data_buffers(&full_raw_sender, &empty_raw_receiver, &mut data_reader) {
                break;
            }
        }
        decoder_thread.join().unwrap();
    })
}
