#![doc(html_favicon_url = "https://example.com/favicon.ico")]
#![doc(html_logo_url = "SAGC-logo-hover.png")]
#![doc(issue_tracker_base_url = "https://github.com/rust-lang/rust/issues/")]

extern crate argparse;

use std::time::Instant;
use std::collections::HashMap;
use std::str;
use chrono;
use argparse::{ArgumentParser, Store, StoreTrue, ParseList};
use mgikit::*;
use termion::terminal_size;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_logo(){
    //let (width, _) = terminal_size().unwrap();
    let (width, _) = match terminal_size() {
        Ok(value) => {
            value
        }
        Err(_) => {
            (80, 0)
        }
    };
    //let width = 80;
    if width > 43{
        //println!("Terminal size: {} columns x {} rows", width, height);
        println! ("\n{}███╗░░░███╗░██████╗░██╗██╗░░██╗██╗████████╗", " ".repeat((width - 43) as usize / 2));
        println! ("{}████╗░████║██╔════╝░██║██║░██╔╝██║╚══██╔══╝", " ".repeat((width - 43) as usize / 2));
        println! ("{}██╔████╔██║██║░░██╗░██║█████═╝░██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! ("{}██║╚██╔╝██║██║░░╚██╗██║██╔═██╗░██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! ("{}██║░╚═╝░██║╚██████╔╝██║██║░╚██╗██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! ("{}╚═╝░░░░░╚═╝░╚═════╝░╚═╝╚═╝░░╚═╝╚═╝░░░╚═╝░░░", " ".repeat((width - 43) as usize / 2));
        println! ("{}                (v{})                 ", " ".repeat((width - 43) as usize / 2), VERSION);
        println! ("{}Ziad Al Bkhetan (ziad.albkhetan@gmail.com)", " ".repeat((width - 43) as usize / 2));
        println! ("{}ACGCGAGACGAGAGATTACGAGCGAGCGAGAGAGACGCTCGAA\n\n\n", " ".repeat((width - 43) as usize / 2));        
    }else{
        println! ("\nMGIKIT\n\n");
    }
}


fn main() {

    print_logo();

    let start = Instant::now();
    println!("Exection start time: {:?}", chrono::offset::Local::now());
    let mut command = String::new();
    let mut arg_read2_file_path = String::new();
    let mut arg_input_folder_path = String::new(); 
    let mut arg_read1_file_path = String::new();
    let mut arg_sample_sheet_file_path = String::new();
    let mut arg_ouput_dir = String::new();

    let mut arg_report_dir = String::new();
    let mut arg_allowed_mismatches = 1;
    let mut arg_template = String::new();
    let mut arg_i7_rc = false;
    let mut arg_i5_rc = false;

    let mut arg_lane = String::new();

    let mut arg_instrument = String::new(); 
    let mut arg_run = String::new(); 
    let mut arg_disable_illumina_format = false;
    let mut arg_keep_barcode = false;

    let mut arg_writing_threshold = 1000;
    let mut arg_read_merging_threshold = 10000;
    //let mut arg_extract_umi = false;

    let mut arg_comprehensive_scan = false;
    let mut arg_undetermined_label = String::from("Undetermined");
    let mut arg_ambiguous_label = String::from("Ambiguous");

    let mut help_messages: HashMap<&str, &str> = HashMap::new();
    let mut arg_force = false;
    let mut arg_report_limit: usize = 50;

    let mut arg_barcode_length: usize = 0;
    let mut arg_testing_reads: usize = 5000;
    let mut arg_popular_template: bool = false;
    let mut arg_no_umi: bool = false;
    let mut arg_max_umi_length: usize = 10;
    
    //QC reports
    let mut arg_qc_report_path : Vec<String> = Vec::new();
    

    help_messages.insert("out_folder", "Path to the output folder. If not provided, the output will be written at mgiKit_ followed by current data and time.");
    help_messages.insert("report_folder", "Prefix of report file. If not provided, the output will be written at output_ followed by current data and time.");
    help_messages.insert(
        "r1_input_path",
        "The path to read1.fastq.gz See the example for the required format.",
    );
    help_messages.insert(
        "r2_input_path",
        "The path to read2.fastq.gz See the example for the required format.",
    );
    help_messages.insert("index_sample_map_file", "The path to the sample/index map.");
    help_messages.insert(
        "template",
        "The general template of the indexes to be used for demultiplexing.",
    );
    help_messages.insert(
        "allowed_mismatches",
        "The number of allowed mismatches when detecting indexes from reads",
    );
    help_messages.insert(
        "demultiplexing",
        "This command is to perform demultiplexing on the dataset",
    );
    help_messages.insert(
        "keep_barcode",
        "Keep the barcode at the tail of read sequence. Default is false.",
    );
    help_messages.insert(
        "instrument_id",
        "The Id of the instrument required for Illumina format.",
    );
    help_messages.insert(
        "run_number",
        "The run number, required for Illumina format.",
    );
    help_messages.insert(
        "lane_number",
        "The lane number, required for Illumina format.",
    );
    help_messages.insert(
        "i5_length",
        "The length of I5, required with the command generateTemplate",
    );
    help_messages.insert(
        "i7_length",
        "The length of I7, required with the command generateTemplate",
    );
    help_messages.insert(
        "testing_reads",
        "Number of reads to perform generateTemplate command, default is 100,000.",
    );
    help_messages.insert("illumina_formate_dis", "Disable illumina file naming and read header format. Output file names and reads' header using MGI format.");
    help_messages.insert(
        "reverse_complement_i7",
        "Convert i7 to reveres complement. Only valid when using general template.",
    );
    help_messages.insert(
        "reverse_complement_i5",
        "Convert i5 to reveres complement. Only valid when using general template.",
    );
    //help_messages.insert("extract_umi", "Extract UMI from reads and attach it the header");
    help_messages.insert("single_read", "Single read input. No read1 file.");
    help_messages.insert(
        "ambiguous_label",
        "The name of the file that contains ambiguous reads.",
    );
    help_messages.insert(
        "undetermined_label",
        "The name of the file that contains undetermined reads.",
    );
    help_messages.insert(
        "merging_threshold",
        "The number of reads that will be merged in one string before writng.",
    );
    help_messages.insert(
        "writing_threshold",
        "The number of merged reads that when reached, data will be saved.",
    );
    help_messages.insert(
        "comprehensive_scan",
        "Check all possible matches, otherwise, search will stop after first match.",
    );
    help_messages.insert(
        "force",
        "Force running th etoll and overwrite existed output/report files.",
    );
    help_messages.insert(
        "report_limit",
        "The number of barcodes to be reported in the list of undetermined and ambiguous barcodes for short/multiqc report. 50 barcodes is the default.",
    );
    help_messages.insert("barcode_length", "The barcode length to detect the template. Default is the length difference between R2 and R1.");
    
    let testing_reads_msg: String = format!("The number of reads used to detect the barcode. Default is {}", arg_testing_reads);
    
    help_messages.insert("testing_reads", &testing_reads_msg);
    help_messages.insert("no_umi", "Don't extract UMI from the read barcode. Default is false.");
    help_messages.insert("popular_template", "Use the most frequest template for all samples even if some of them have more matches with other template. Default is true.");
    help_messages.insert("max_umi_length", "The maximum expected UMI length. Default is 10.");
    
    help_messages.insert("", "");

    {
        // this block limits scope of borrows by ap.refer() method
        let mut ap = ArgumentParser::new();
        ap.set_description("MGIKit - Demultiplexer.");

        ap.refer(&mut command).add_argument(
            "command",
            Store,
            help_messages.get("demultiplexing").unwrap(),
        );
        
       

        //ap.add_option(&["-V", "--version"],
        //    println!(env!("CARGO_PKG_VERSION").to_string()),
        //    help_messages.get("").unwrap());
        
        ap.refer(&mut arg_input_folder_path).add_option(
            &["-i", "--input"],
            Store,
            help_messages.get("r2_input_path").unwrap(),
        );

        ap.refer(&mut arg_qc_report_path).add_option(
            &["--qc-reports"],
            ParseList,
            help_messages.get("r2_input_path").unwrap(),
        );

        ap.refer(&mut arg_read2_file_path).add_option(
            &["-r", "--read2", "-2"],
            Store,
            help_messages.get("r2_input_path").unwrap(),
        );

        ap.refer(&mut arg_read1_file_path).add_option(
            &["-f", "--read1", "-1"],
            Store,
            help_messages.get("r1_input_path").unwrap(),
        );
        ap.refer(&mut arg_sample_sheet_file_path).add_option(
            &["-s", "--sample-sheet"],
            Store,
            help_messages.get("index_sample_map_file").unwrap(),
        );

        ap.refer(&mut arg_ouput_dir).add_option(
            &["-o", "--output"],
            Store,
            help_messages.get("out_folder").unwrap(),
        );

        ap.refer(&mut arg_report_dir).add_option(
            &["--reports"],
            Store,
            help_messages.get("report_folder").unwrap(),
        );

        ap.refer(&mut arg_allowed_mismatches).add_option(
            &["-m", "--mismatches"],
            Store,
            help_messages.get("allowed_mismatches").unwrap(),
        );

        ap.refer(&mut arg_template).add_option(
            &["--template"],
            Store,
            help_messages.get("template").unwrap(),
        );

        ap.refer(&mut arg_i7_rc).add_option(
            &["--i7-rc"],
            StoreTrue,
            help_messages.get("reverse_complement_i7").unwrap(),
        );

        ap.refer(&mut arg_i5_rc).add_option(
            &["--i5-rc"],
            StoreTrue,
            help_messages.get("reverse_complement_i5").unwrap(),
        );

        ap.refer(&mut arg_disable_illumina_format).add_option(
            &["--disable-illumina"],
            StoreTrue,
            help_messages.get("illumina_formate_dis").unwrap(),
        );

        ap.refer(&mut arg_keep_barcode).add_option(
            &["--keep-barcode"],
            StoreTrue,
            help_messages.get("keep_barcode").unwrap(),
        );

        ap.refer(&mut arg_writing_threshold).add_option(
            &["--writing-buffer"],
            Store,
            help_messages.get("writing_threshold").unwrap(),
        );

        ap.refer(&mut arg_read_merging_threshold).add_option(
            &["--merged-reads"],
            Store,
            help_messages.get("merging_threshold").unwrap(),
        );

        ap.refer(&mut arg_lane).add_option(
            &["--lane"],
            Store,
            help_messages.get("lane_number").unwrap(),
        );

        ap.refer(&mut arg_instrument).add_option(
            &["--instrument"],
            Store,
            help_messages.get("instrument_id").unwrap(),
        );

        ap.refer(&mut arg_run).add_option(
            &["--run"],
            Store,
            help_messages.get("run_number").unwrap(),
        );

        ap.refer(&mut arg_undetermined_label).add_option(
            &["--undetermined-label"],
            Store,
            help_messages.get("undetermined_label").unwrap(),
        );

        ap.refer(&mut arg_ambiguous_label).add_option(
            &["--ambiguous-label"],
            Store,
            help_messages.get("ambiguous_label").unwrap(),
        );

        ap.refer(&mut arg_comprehensive_scan).add_option(
            &["--comprehensive-scan"],
            StoreTrue,
            help_messages.get("comprehensive_scan").unwrap(),
        );

        ap.refer(&mut arg_force).add_option(
            &["--force"],
            StoreTrue,
            help_messages.get("force").unwrap(),
        );
        ap.refer(&mut arg_report_limit).add_option(
            &["--report-limit"],
            Store,
            help_messages.get("report_limit").unwrap(),
        );
        ap.refer(&mut arg_barcode_length).add_option(
            &["--barcode-length"],
            Store,
            help_messages.get("barcode_length").unwrap(),
        );


        
    
    
        ap.refer(&mut arg_popular_template).add_option(
            &["--popular-template"],
            StoreTrue,
            help_messages.get("popular_template").unwrap(),
        );

        ap.refer(&mut arg_no_umi).add_option(
            &["--no-umi"],
            StoreTrue,
            help_messages.get("no_umi").unwrap(),
        );
            
        ap.refer(&mut arg_testing_reads).add_option(
            &["--testing-reads"],
            Store,
            help_messages.get("testing_reads").unwrap(),
        );
        
        ap.refer(&mut arg_max_umi_length).add_option(
            &["--max-umi-len"],
            Store,
            help_messages.get("max_umi_length").unwrap(),
        );

        ap.parse_args_or_exit();  
    }

    if command == "demultiplex" {
        println!("Perform demutiplex command");

        demultiplex(
            &arg_input_folder_path,
            &mut arg_read1_file_path,
            &mut arg_read2_file_path,
            &arg_sample_sheet_file_path,
            &arg_ouput_dir,
            &arg_report_dir,
            arg_allowed_mismatches,
            &arg_template,
            arg_i7_rc,
            arg_i5_rc,
            &arg_lane,
            &arg_instrument,
            &arg_run,
            arg_disable_illumina_format,
            arg_keep_barcode,
            arg_writing_threshold,
            arg_read_merging_threshold,
            arg_comprehensive_scan,
            &arg_undetermined_label,
            &arg_ambiguous_label,
            arg_force,
            arg_report_limit
        );
    } else if command == "reports" {
        println!("Perform reports command");
        merge_qc_reports(&arg_qc_report_path, &arg_ouput_dir);
    }else if command == "template" {
        println!("Perform detect template command");
        detect_template(&arg_read1_file_path, 
            &arg_read2_file_path,
            &arg_sample_sheet_file_path,
            &arg_ouput_dir,
            arg_testing_reads,
            arg_barcode_length,
            ! arg_no_umi,
            arg_popular_template,
            arg_max_umi_length
        );
    }else {
        panic!("Please enter a command to perform from (demultiplex, reports, and template)!");
    }

    let dur = start.elapsed();
    println!("{} seconds for performing the task.", dur.as_secs());
    println!("Exection end time: {:?}", chrono::offset::Local::now());
    //Ok()
}

