use std::fs;
use std::path::{ Path, PathBuf };
use log::info;
use core::panic;
use flate2::read::MultiGzDecoder;
use std::io::{ BufReader, BufWriter, Read, Write };
use std::fs::File;
use std::thread::{ self, JoinHandle };
use crossbeam_channel::{ Receiver, Sender };
use memchr::memchr_iter;
use log::debug;

fn read_bytes_in_reads<R: Read>(
    reader: &mut R,
    buffer: &mut [u8],
    _minimum: usize,
    last_byte: &mut usize
) -> (bool, Vec<usize>) {
    let mut curr_bytes: usize;
    loop {
        //debug!("buffer length: {}, starting from {}", buffer.len(), last_byte);
        curr_bytes = reader.read(&mut buffer[*last_byte..]).unwrap();
        if curr_bytes == 0 {
            //debug!("total lines: {} - no more", line_cnt);
            break;
        }
        *last_byte += curr_bytes;
    }
    //debug!("total lines: {} - still more", line_cnt);
    return (*last_byte > 0, memchr_iter(b'\n', &buffer[..*last_byte]).collect::<Vec<usize>>());
}

pub fn fill_send_buffers<R: Read>(
    full_sender: &Sender<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver: &Receiver<(usize, Vec<u8>)>,
    reader: &mut R,
    read_cnt: usize,
    extra: &mut Vec<u8>,
    extra_len: &mut usize
) -> bool {
    let mut total_bytes = 0;
    match empty_receiver.recv() {
        Ok((_, mut buffer)) => {
            if *extra_len > 0 {
                buffer[..*extra_len].copy_from_slice(&extra[..*extra_len]);
                total_bytes = *extra_len;
                *extra_len = 0;
            }

            let (_, mut lines) = read_bytes_in_reads(
                reader,
                &mut buffer[..],
                read_cnt,
                &mut total_bytes
            );
            //let mut lines = memchr_iter(b'\n', &buffer[..total_bytes]);
            //debug!("total read bytes: {}  -  lines: {}  - needed lines: {}", total_bytes, lines.count(), read_cnt * 4);
            //let mut lines = memchr_iter(b'\n', &buffer[..total_bytes]);
            if lines.len() > read_cnt * 4 {
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

pub fn parallel_reader_thread(
    paired_reads: PathBuf,
    barcode_reads: PathBuf,
    batch_size: usize,
    read_rb: bool,
    read_rp: bool,
    full_sender_rb: Sender<(usize, Vec<u8>, Vec<usize>)>,
    full_sender_rp: Sender<(usize, Vec<u8>, Vec<usize>)>,
    full_receiver_rb: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    full_receiver_rp: Receiver<(usize, Vec<u8>, Vec<usize>)>,
    empty_receiver_rb: Receiver<(usize, Vec<u8>)>,
    empty_receiver_rp: Receiver<(usize, Vec<u8>)>,
    full_sender: Sender<(usize, Vec<u8>, Vec<usize>, usize, Vec<u8>, Vec<usize>)>,
    processing_threads: usize,
    paired_input: bool
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

        let mut extra_rb: Vec<u8> = if read_rb { vec![b'0'; 1000000] } else { Vec::new() };
        let mut extra_len_rb: usize = 0;

        let mut extra_rp: Vec<u8> = if read_rp { vec![b'0'; 1000000] } else { Vec::new() };
        let mut extra_len_rp: usize = 0;

        let mut keep_reading_rp = read_rp;
        let mut readers_finished = false;
        loop {
            if read_rb {
                match reader_barcode_read {
                    Some(ref mut reader_barcode) => {
                        let _ = fill_send_buffers(
                            &full_sender_rb,
                            &empty_receiver_rb,
                            reader_barcode,
                            batch_size,
                            &mut extra_rb,
                            &mut extra_len_rb
                        );
                        debug!("Sent rb full buffer!");
                    }
                    None => {}
                }
            }

            if read_rp {
                match reader_paired_read {
                    Some(ref mut reader_paired) => {
                        keep_reading_rp = fill_send_buffers(
                            &full_sender_rp,
                            &empty_receiver_rp,
                            reader_paired,
                            batch_size,
                            &mut extra_rp,
                            &mut extra_len_rp
                        );
                        debug!("Sent rp full buffer!");
                    }
                    None => {}
                }
                if read_rp && !read_rb && !keep_reading_rp {
                    readers_finished = true;
                }
            }

            if read_rb {
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
                if read_rb {
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
