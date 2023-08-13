#![doc(html_favicon_url = "https://example.com/favicon.ico")]
#![doc(html_logo_url = "assets/SAGC-logo-hover.png")]
#![doc(issue_tracker_base_url = "https://github.com/sagc-bioinformatics/mgikit/issues")]

use std::time::Instant;
use std::str;
use chrono;
use mgikit::*;
use termion::terminal_size;
use clap::{ArgAction, Command, Arg};

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
        println! (  "{}████╗░████║██╔════╝░██║██║░██╔╝██║╚══██╔══╝", " ".repeat((width - 43) as usize / 2));
        println! (  "{}██╔████╔██║██║░░██╗░██║█████═╝░██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! (  "{}██║╚██╔╝██║██║░░╚██╗██║██╔═██╗░██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! (  "{}██║░╚═╝░██║╚██████╔╝██║██║░╚██╗██║░░░██║░░░", " ".repeat((width - 43) as usize / 2));
        println! (  "{}╚═╝░░░░░╚═╝░╚═════╝░╚═╝╚═╝░░╚═╝╚═╝░░░╚═╝░░░", " ".repeat((width - 43) as usize / 2));
        println! (  "{}                (v{})                 ", " ".repeat((width - 43) as usize / 2), VERSION);
        println! (  "{}Ziad Al-Bkhetan (ziad.albkhetan@gmail.com)", " ".repeat((width - 43) as usize / 2));
        println! (  "{}ACGCGAGACGAGAGATTACGAGCGAGCGAGAGAGACGCTCGAA\n\n\n", " ".repeat((width - 43) as usize / 2));        
    }else{
        println! ("\nMGIKIT\n\n");
    }
}


fn main() {

    print_logo();

    {
        let matches = Command::new("MGIKIT - MGI data demultipexing kit.")
        .about("mgikit is a multiple commands to support MGI fastq demultiplexing.")
        .author("Ziad Al-Bkhetan, ziad.albkhetan@gmail.com")
        .version(VERSION)
        .propagate_version(true)
        .after_help("Check out mgikit documenation at `https://github.com/sagc-bioinformatics/mgikit` for more details.")
        .subcommand(
            Command::new("demultiplex")
                .about("Demultipex fastq files.")
                .arg(
                    Arg::new("arg_input_folder_path")
                        .short('i')
                        .long("input")
                        .help("The path to read2.fastq.gz See the example for the required format.")
                )
                .arg(
                    Arg::new("arg_read2_file_path")
                        .short('r')
                        .long("read2")
                        .alias("2")
                        .help("The path to read2.fastq.gz See the example for the required format.")
                )
                .arg(
                    Arg::new("arg_read1_file_path")
                        .short('f')
                        .long("read1")
                        .alias("1")
                        .help("The path to read1.fastq.gz See the example for the required format.")
                )
                .arg(
                    Arg::new("arg_sample_sheet_file_path")
                        .short('s')
                        .long("sample-sheet")
                        .required(true)
                        .help("The path to the sample/index map.")
                )
                .arg(
                    Arg::new("arg_ouput_dir")
                        .short('o')
                        .long("output")
                        .help("Path to the output folder. If not provided, the output will be written at mgiKit_ followed by current data and time.")
                )
                .arg(
                    Arg::new("arg_report_dir")
                        .long("reports")
                        .help("Prefix of report file. If not provided, the output will be written at output_ followed by current data and time.")
                )
                .arg(
                    Arg::new("arg_allowed_mismatches")
                        .short('m')
                        .long("mismatches")
                        .default_value("1")
                        .value_parser(clap::value_parser!(usize))
                        .help("The number of allowed mismatches when detecting indexes from reads.")
                )
                .arg(
                    Arg::new("arg_template")
                        .long("template")
                        
                        .help("The general template of the indexes to be used for demultiplexing.")
                )
                .arg(
                    Arg::new("arg_i7_rc")
                        .long("i7-rc")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Convert i7 to reveres complement. Only valid when using general template.")
                )
                .arg(
                    Arg::new("arg_i5_rc")
                        .long("i5-rc")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Convert i5 to reveres complement. Only valid when using general template.")
                )
                .arg(
                    Arg::new("arg_disable_illumina_format")
                        .long("disable-illumina")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Disable illumina file naming and read header format. Output file names and reads' header using MGI format.")
                )
                .arg(
                    Arg::new("arg_keep_barcode")
                        .long("keep-barcode")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Keep the barcode at the tail of read sequence.")
                )
                .arg(
                    Arg::new("arg_writing_threshold")
                        .long("writing-buffer")
                        .default_value("10000")
                        .value_parser(clap::value_parser!(usize))
                        .help("The number of merged reads that when reached, data will be saved.")
                )
                .arg(
                    Arg::new("arg_read_merging_threshold")
                        .long("merged-reads")
                        .default_value("1000")
                        .value_parser(clap::value_parser!(usize))
                        .help("The number of reads that will be merged in one string before writng.")
                )
                .arg(
                    Arg::new("arg_lane")
                        .long("lane")
                        
                        .help("The lane number, required for Illumina format.")
                )
                .arg(
                    Arg::new("arg_instrument")
                        .long("instrument")
                        
                        .help("The Id of the instrument required for Illumina format.")
                )
                .arg(
                    Arg::new("arg_run")
                        .long("run")
                        
                        .help("The run number, required for Illumina format.")
                )
                .arg(
                    Arg::new("arg_undetermined_label")
                        .long("undetermined-label")
                        .default_value("Undetermined")
                        .help("The name of the file that contains undetermined reads.")
                )
                .arg(
                    Arg::new("arg_ambiguous_label")
                        .long("ambiguous-label")
                        .default_value("Ambiguous")
                        .help("The name of the file that contains ambiguous reads.")
                )
                .arg(
                    Arg::new("arg_comprehensive_scan")
                        .long("comprehensive-scan")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Check all possible matches, otherwise, search will stop after first match. Only needed for mixed library.")
                )
                .arg(
                    Arg::new("arg_force")
                        .long("force")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Force running th etoll and overwrite existed output/report files.")
                )
                .arg(
                    Arg::new("arg_report_limit")
                        .long("report-limit")
                        .default_value("20")
                        .value_parser(clap::value_parser!(usize))
                        .help("The number of barcodes to be reported in the list of undetermined and ambiguous barcodes for short/multiqc report.")
                )
                .arg(
                    Arg::new("arg_read1_file_name_suf")
                        .long("r1-file-suf")
                        .default_value("_read_1.fq.gz")
                        .help("The suffix to read1 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read1 file. There should be one file with this suffix in the input directory.")
                )
                .arg(
                    Arg::new("arg_read2_file_name_suf")
                        .long("r2-file-suf")
                        .default_value("_read_2.fq.gz")
                        .help("The suffix to read2 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read2 file. There should be one file with this suffix in the input directory.")
                )
                .arg(
                    Arg::new("arg_info_file")
                        .long("info-file")
                        .default_value("BioInfo.csv")
                        .help("The name of the info file that contains the run information. Only needed when using the `--input` parameter.")
                )
        )
        .subcommand(
            Command::new("template")
                .about("Detect barcode template.")
                .arg(
                    Arg::new("arg_barcode_length")
                        .long("barcode-length")
                        .default_value("0")
                        .value_parser(clap::value_parser!(usize))
                        .help("The barcode length to detect the template. When set to 0, the barcode length will be length's difference between R2 and R1.")
                )
                .arg(
                    Arg::new("arg_popular_template")
                        .long("popular-template")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Use the most frequest template for all samples even if some of them have more matches with other template.")
                )
                .arg(
                    Arg::new("arg_no_umi")
                        .long("no-umi")
                        .action(ArgAction::SetTrue)
                        .default_value("false")
                        .help("Don't extract UMI from the read barcode.")
                )
                .arg(
                    Arg::new("arg_testing_reads")
                        .long("testing-reads")
                        .default_value("5000")
                        .value_parser(clap::value_parser!(usize))
                        .help("The number of reads used to detect the barcode.")
                )
                .arg(
                    Arg::new("arg_max_umi_length")
                        .long("max-umi-len")
                        .default_value("10")
                        .value_parser(clap::value_parser!(usize))
                        .help("The maximum expected UMI length.")
                )
        )
        .subcommand(
            Command::new("report")
                .about("Merge demultipexing reports.")
                .arg(
                    Arg::new("arg_qc_report_path")
                        .long("qc-reports")
                        .required(true)
                        .value_delimiter(' ') 
                        .help("The paths to the QC reports separated by a space ` `.")
                ) 
                .arg(
                    Arg::new("arg_ouput_dir")
                        .short('o')
                        .long("output")
                        .help("output directory")
                )               
        )
        .get_matches();
        
        let start = Instant::now();
        println!("Exection start time: {:?}", chrono::offset::Local::now());
    
        match matches.subcommand() {
            Some(("demultiplex", demultiplex_command)) => {
                let arg_input_folder_path: &String = demultiplex_command.get_one::<String>("arg_input_folder_path").unwrap();
                let mut arg_read1_file_path: String = demultiplex_command.get_one::<String>("arg_read1_file_path").unwrap().to_string();
                let mut arg_read2_file_path: String = demultiplex_command.get_one::<String>("arg_read2_file_path").unwrap().to_string();
                let arg_sample_sheet_file_path: &String = demultiplex_command.get_one::<String>("arg_sample_sheet_file_path").unwrap();
                let arg_ouput_dir: &String = demultiplex_command.get_one::<String>("arg_ouput_dir").unwrap();
                let arg_report_dir: &String = demultiplex_command.get_one::<String>("arg_report_dir").unwrap();
                let arg_allowed_mismatches: &usize = demultiplex_command.get_one::<usize>("arg_allowed_mismatches").unwrap();
                let arg_template: &String = demultiplex_command.get_one::<String>("arg_template").unwrap();
                let arg_i7_rc: &bool = demultiplex_command.get_one::<bool>("arg_i7_rc").unwrap();
                let arg_i5_rc: &bool = demultiplex_command.get_one::<bool>("arg_i5_rc").unwrap();
                let arg_lane: &String = demultiplex_command.get_one::<String>("arg_lane").unwrap();
                let arg_instrument: &String = demultiplex_command.get_one::<String>("arg_instrument").unwrap();
                let arg_run: &String = demultiplex_command.get_one::<String>("arg_run").unwrap();
                let arg_disable_illumina_format: &bool = demultiplex_command.get_one::<bool>("arg_disable_illumina_format").unwrap();
                let arg_keep_barcode: &bool = demultiplex_command.get_one::<bool>("arg_keep_barcode").unwrap();
                let arg_writing_threshold: &usize = demultiplex_command.get_one::<usize>("arg_writing_threshold").unwrap();
                let arg_read_merging_threshold: &usize = demultiplex_command.get_one::<usize>("arg_read_merging_threshold").unwrap();
                let arg_comprehensive_scan: &bool = demultiplex_command.get_one::<bool>("arg_comprehensive_scan").unwrap();
                let arg_undetermined_label: &String = demultiplex_command.get_one::<String>("arg_undetermined_label").unwrap();
                let arg_ambiguous_label: &String = demultiplex_command.get_one::<String>("arg_ambiguous_label").unwrap();
                let arg_force: &bool = demultiplex_command.get_one::<bool>("arg_force").unwrap();
                let arg_report_limit: &usize = demultiplex_command.get_one::<usize>("arg_report_limit").unwrap();
                let arg_read1_file_name_suf: &String = demultiplex_command.get_one::<String>("arg_read1_file_name_suf").unwrap();
                let arg_read2_file_name_suf: &String = demultiplex_command.get_one::<String>("arg_read2_file_name_suf").unwrap();
                let arg_info_file: &String = demultiplex_command.get_one::<String>("arg_info_file").unwrap();
                
                demultiplex(
                    arg_input_folder_path,
                    &mut arg_read1_file_path,
                    &mut arg_read2_file_path,
                    arg_sample_sheet_file_path,
                    arg_ouput_dir,
                    arg_report_dir,
                    *arg_allowed_mismatches,
                    arg_template,
                    *arg_i7_rc,
                    *arg_i5_rc,
                    arg_lane,
                    arg_instrument,
                    arg_run,
                    *arg_disable_illumina_format,
                    *arg_keep_barcode,
                    *arg_writing_threshold,
                    *arg_read_merging_threshold,
                    *arg_comprehensive_scan,
                    arg_undetermined_label,
                    arg_ambiguous_label,
                    *arg_force,
                    *arg_report_limit,
                    arg_read1_file_name_suf,
                    arg_read2_file_name_suf,
                    arg_info_file
                );
            },
            Some(("report", report_command)) => {
                let arg_ouput_dir: &String = report_command.get_one::<String>("arg_ouput_dir").unwrap();
                let arg_qc_report_path: Vec<String> = report_command.get_many::<String>("arg_qc_report_path").unwrap().map(|it: &String| it.to_string()).collect::<Vec<_>>();
                merge_qc_reports(
                    &arg_qc_report_path, 
                    arg_ouput_dir);
                
            },
            Some(("template", template_command)) => {
                let arg_read1_file_path: String = template_command.get_one::<String>("arg_read1_file_path").unwrap().to_string();
                let arg_read2_file_path: String = template_command.get_one::<String>("arg_read2_file_path").unwrap().to_string();
                let arg_sample_sheet_file_path: &String = template_command.get_one::<String>("arg_sample_sheet_file_path").unwrap();
                let arg_barcode_length: &usize = template_command.get_one::<usize>("arg_barcode_length").unwrap();
                let arg_testing_reads: &usize = template_command.get_one::<usize>("arg_testing_reads").unwrap();
                let arg_max_umi_length: &usize = template_command.get_one::<usize>("arg_max_umi_length").unwrap();
                let arg_popular_template: &bool = template_command.get_one::<bool>("arg_popular_template").unwrap();
                let arg_no_umi: &bool = template_command.get_one::<bool>("arg_no_umi").unwrap();
                let arg_ouput_dir: &String = template_command.get_one::<String>("arg_ouput_dir").unwrap();
                
                detect_template(
                    &arg_read1_file_path, 
                    &arg_read2_file_path,
                    arg_sample_sheet_file_path,
                    arg_ouput_dir,
                    *arg_testing_reads,
                    *arg_barcode_length,
                    ! arg_no_umi,
                    *arg_popular_template,
                    *arg_max_umi_length
                );
            },
            Some((command_nm, _)) => {
                panic!("Unknown command `{}`. Please enter a command to perform from (demultiplex, report, or template)!", command_nm);
            }
            None => {
                panic!("Please enter a command to perform from (demultiplex, report, or template)!");
            }
        }
        let dur = start.elapsed();
        println!("{} seconds for performing the task.", dur.as_secs());
        println!("Exection end time: {:?}", chrono::offset::Local::now());
    }
}

