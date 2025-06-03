#![doc(html_favicon_url = "https://example.com/favicon.ico")]
#![doc(html_logo_url = "assets/SAGC-logo-hover.png")]
#![doc(issue_tracker_base_url = "https://github.com/sagc-bioinformatics/mgikit/issues")]

use std::time::Instant;
use std::str;
use chrono;
use mgikit::*;
use termion::terminal_size;
use clap::{ ArgAction, Command, Arg };
use log::{ warn, info, error, LevelFilter };
use std::env;
use env_logger::{ Builder, Target };
use std::str::FromStr;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_logo() {
    //let (width, _) = terminal_size().unwrap();
    let (width, _) = match terminal_size() {
        Ok(value) => { value }
        Err(_) => { (80, 0) }
    };
    //let width = 80;
    if width > 43 {
        //println!("Terminal size: {} columns x {} rows", width, height);
        println!(
            "\n{}███╗░░░███╗░██████╗░██╗██╗░░██╗██╗████████╗",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}████╗░████║██╔════╝░██║██║░██╔╝██║╚══██╔══╝",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}██╔████╔██║██║░░██╗░██║█████═╝░██║░░░██║░░░",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}██║╚██╔╝██║██║░░╚██╗██║██╔═██╗░██║░░░██║░░░",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}██║░╚═╝░██║╚██████╔╝██║██║░╚██╗██║░░░██║░░░",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}╚═╝░░░░░╚═╝░╚═════╝░╚═╝╚═╝░░╚═╝╚═╝░░░╚═╝░░░",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}                (v{})                 ",
            " ".repeat(((width - 43) as usize) / 2),
            VERSION
        );
        println!(
            "{}Ziad Al-Bkhetan (ziad.albkhetan@gmail.com)",
            " ".repeat(((width - 43) as usize) / 2)
        );
        println!(
            "{}ACGCGAGACGAGAGATT-MGIKIT-GCGAGAGAGACGCTCGAA\n\n\n",
            " ".repeat(((width - 43) as usize) / 2)
        );
    } else {
        println!("\nMGIKIT\n\n");
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
            .after_help(
                "Check out mgikit documenation at `https://github.com/sagc-bioinformatics/mgikit` for more details."
            )
            .arg(
                Arg::new("arg_log_level")
                    .long("log-level")
                    .short('l')
                    .global(true)
                    .help("Set the log level (error, warn, info, debug, trace)")
                    .default_value("info")
                    .help(
                        "log level for output messages. Expected values: [error, warn, info, debug, trace]. Default is info."
                    )
            )
            .subcommand(
                Command::new("demultiplex")
                    .about("Demultipex fastq files.")
                    .arg(
                        Arg::new("arg_input_folder_path")
                            .short('i')
                            .default_value("")
                            .long("input")
                            .help(
                                "The path to read2.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_read2_file_path")
                            .short('r')
                            .long("read2")
                            .alias("2")
                            .default_value("")
                            .help(
                                "The path to read2.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_read1_file_path")
                            .short('f')
                            .long("read1")
                            .alias("1")
                            .default_value("")
                            .help(
                                "The path to read1.fastq.gz See the example for the required format."
                            )
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
                            .default_value("")
                            .help(
                                "Path to the output folder. If not provided, the output will be written at mgiKit_ followed by current data and time."
                            )
                    )
                    .arg(
                        Arg::new("arg_report_dir")
                            .long("reports")
                            .default_value("")
                            .help(
                                "Prefix of report file. If not provided, the output will be written at output_ followed by current data and time."
                            )
                    )
                    .arg(
                        Arg::new("arg_allowed_mismatches")
                            .short('m')
                            .long("mismatches")
                            .default_value("1")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The number of allowed mismatches when detecting indexes from reads."
                            )
                    )
                    .arg(
                        Arg::new("arg_template")
                            .long("template")
                            .default_value("")
                            .help(
                                "The general template of the indexes to be used for demultiplexing."
                            )
                    )
                    .arg(
                        Arg::new("arg_i7_rc")
                            .long("i7-rc")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Convert i7 to reveres complement. Only valid when using general template."
                            )
                    )
                    .arg(
                        Arg::new("arg_i5_rc")
                            .long("i5-rc")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Convert i5 to reveres complement. Only valid when using general template."
                            )
                    )
                    .arg(
                        Arg::new("arg_disable_illumina_format")
                            .long("disable-illumina")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Disable illumina file naming and read header format. Output file names and reads' header using MGI format."
                            )
                    )
                    .arg(
                        Arg::new("arg_keep_barcode")
                            .long("keep-barcode")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help("Keep the barcode at the tail of read sequence.")
                    )
                    .arg(
                        Arg::new("arg_writing_buffer_size")
                            .long("writing-buffer-size")
                            .default_value("67108864")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The size of the buffer for each sample to be filled with data then written once to the disk."
                            )
                    )
                    .arg(
                        Arg::new("arg_lane")
                            .long("lane")
                            .default_value("")
                            .help("The lane number, required for Illumina format.")
                    )
                    .arg(
                        Arg::new("arg_instrument")
                            .long("instrument")
                            .default_value("")
                            .help("The Id of the instrument required for Illumina format.")
                    )
                    .arg(
                        Arg::new("arg_run")
                            .long("run")
                            .default_value("")
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
                            .help(
                                "Check all possible matches, otherwise, search will stop after first match. Only needed for mixed library."
                            )
                    )
                    .arg(
                        Arg::new("arg_force")
                            .long("force")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Force running the tool and overwrite existed output/report files."
                            )
                    )
                    .arg(
                        Arg::new("arg_report_limit")
                            .long("report-limit")
                            .default_value("20")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The number of barcodes to be reported in the list of undetermined and ambiguous barcodes for short/multiqc report."
                            )
                    )
                    .arg(
                        Arg::new("arg_read1_file_name_suf")
                            .long("in-r1-file-suf")
                            .default_value("_read_1.fq.gz")
                            .help(
                                "The suffix to read1 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read1 file. There should be one file with this suffix in the input directory."
                            )
                    )
                    .arg(
                        Arg::new("arg_read2_file_name_suf")
                            .long("in-r2-file-suf")
                            .default_value("_read_2.fq.gz")
                            .help(
                                "The suffix to read2 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read2 file. There should be one file with this suffix in the input directory."
                            )
                    )
                    .arg(
                        Arg::new("arg_info_file")
                            .long("info-file")
                            .default_value("")
                            .help(
                                "The path to the info file that contains the run information (similar to `BioInfo.csv` generated by MGI machines under the lane directory). Check the documentation for more details."
                            )
                    )
                    .arg(
                        Arg::new("arg_report_level")
                            .long("report-level")
                            .default_value("2")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The level of reporting. 0 no reports will be generated!, 1 data quality and demultipexing reports. 2: all reports (reports on data quality, demultipexing, undetermined and ambigouse barcodes)."
                            )
                    )
                    .arg(
                        Arg::new("arg_compression_level")
                            .long("compression-level")
                            .default_value("1")
                            .value_parser(clap::value_parser!(u32))
                            .help(
                                "The level of compression (between 0 and 12). 0 is fast but no compression, 12 is slow but high compression."
                            )
                    )
                    .arg(
                        Arg::new("arg_compression_buffer_size")
                            .long("compression-buffer-size")
                            .default_value("131072")
                            .value_parser(clap::value_parser!(usize))
                            .help("The size of the buffer for data compression for each sample.")
                    )
                    .arg(
                        Arg::new("arg_ignore_undetermined")
                            .long("ignore-undetermined")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Do not stop if there are many undetermined reads in the dataset."
                            )
                    )
                    .arg(
                        Arg::new("arg_log_path")
                            .long("log")
                            .default_value("")
                            .help("Path to the output log, instead of writing to the stdout.")
                    )
                    .arg(
                        Arg::new("arg_mgi_full_header")
                            .long("mgi-full-header")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Add sample barcode and UMI to MGI header when MGI format is enabled."
                            )
                    )
                    .arg(
                        Arg::new("arg_all_index_error")
                            .long("all-index-error")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "By default, the allowed mismatches `-m or --mismatches` are considered to be per index. This flag will make it for the total mismatches across all indices."
                            )
                    )
                    .arg(
                        Arg::new("arg_memory")
                            .long("memory")
                            .default_value("0")
                            .value_parser(clap::value_parser!(f64))
                            .help(
                                "The requested maximum memory to be used (in giga byte). Check the documentation for memory optimisation options. Default is 0 then the tool will use the available memory on the machine."
                            )
                    )
                    .arg(
                        Arg::new("arg_not_mgi")
                            .long("not-mgi")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "This flag needs to be enabled if the input fastq files don't have MGI format."
                            )
                    )
                    .arg(
                        Arg::new("arg_threads")
                            .long("threads")
                            .short('t')
                            .default_value("0")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The requested threads to be used. Default is 0 which means to use all available CPUs."
                            )
                    )
                    .arg(
                        Arg::new("arg_check_content")
                            .long("validate")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help("Validate the content of the fastq files.")
                    )
            )
            .subcommand(
                Command::new("template")
                    .about("Detect barcode template.")
                    .arg(
                        Arg::new("arg_input_folder_path")
                            .short('i')
                            .default_value("")
                            .long("input")
                            .help(
                                "The path to read2.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_read2_file_path")
                            .short('r')
                            .long("read2")
                            .alias("2")
                            .default_value("")
                            .help(
                                "The path to read2.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_read1_file_path")
                            .short('f')
                            .long("read1")
                            .alias("1")
                            .default_value("")
                            .help(
                                "The path to read1.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_sample_sheet_file_path")
                            .short('s')
                            .long("sample-sheet")
                            .required(true)
                            .help("The path to the sample/index map.")
                    )
                    .arg(
                        Arg::new("arg_info_file")
                            .long("info-file")
                            .default_value("")
                            .help(
                                "The path to the info file that contains the run information (similar to `BioInfo.csv` generated by MGI machines under the lane directory). Check the documentation for more details."
                            )
                    )
                    .arg(
                        Arg::new("arg_ouput_dir")
                            .short('o')
                            .long("output")
                            .default_value("")
                            .help(
                                "Path to the output file prefix. If not provided, the output will be written at mgiKit_ followed by current data and time."
                            )
                    )
                    .arg(
                        Arg::new("arg_barcode_length")
                            .long("barcode-length")
                            .default_value("0")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The barcode length to detect the template. When set to 0, the barcode length will be length's difference between R2 and R1."
                            )
                    )
                    .arg(
                        Arg::new("arg_popular_template")
                            .long("popular-template")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Use the most frequent template for all samples even if some of them have more matches with other template."
                            )
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
                    .arg(
                        Arg::new("arg_read1_file_name_suf")
                            .long("in-r1-file-suf")
                            .default_value("_read_1.fq.gz")
                            .help(
                                "The suffix to read1 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read1 file. There should be one file with this suffix in the input directory."
                            )
                    )
                    .arg(
                        Arg::new("arg_read2_file_name_suf")
                            .long("in-r2-file-suf")
                            .default_value("_read_2.fq.gz")
                            .help(
                                "The suffix to read2 file name. When using the --input parameter, the tool looks for the file that ends with this suffix and use it as read2 file. There should be one file with this suffix in the input directory."
                            )
                    )
            )
            .subcommand(
                Command::new("report")
                    .about("Merge demultipexing reports.")
                    .arg(
                        Arg::new("arg_qc_report_path")
                            .long("qc-report")
                            .required(true)
                            .action(clap::ArgAction::Append)
                            .help("The paths to the QC reports, repeat it for each report.")
                    )
                    .arg(
                        Arg::new("arg_ouput_dir").short('o').long("output").help("output directory")
                    )
                    .arg(
                        Arg::new("arg_lane")
                            .long("lane")
                            .default_value("all")
                            .help("The lane number, required for report name.")
                    )
                    .arg(
                        Arg::new("arg_prefix")
                            .long("prefix")
                            .default_value("")
                            .help(
                                "The prefix of the report. By default, it is the first part of the last input report."
                            )
                    )
            )
            .subcommand(
                Command::new("reformat")
                    .about("Reformat MGI fastq headers to Illumina's and prepare quality report.")
                    .arg(
                        Arg::new("arg_read2_file_path")
                            .short('r')
                            .long("read2")
                            .alias("2")
                            .default_value("")
                            .help(
                                "The path to read2.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_read1_file_path")
                            .short('f')
                            .long("read1")
                            .alias("1")
                            .default_value("")
                            .help(
                                "The path to read1.fastq.gz See the example for the required format."
                            )
                    )
                    .arg(
                        Arg::new("arg_ouput_dir")
                            .short('o')
                            .long("output")
                            .default_value("")
                            .help(
                                "Path to the output folder. If not provided, the output will be written at mgiKit_ followed by current data and time."
                            )
                    )
                    .arg(
                        Arg::new("arg_report_dir")
                            .long("reports")
                            .default_value("")
                            .help(
                                "Prefix of report file. If not provided, the output will be written at output_ followed by current data and time."
                            )
                    )
                    .arg(
                        Arg::new("arg_disable_illumina_format")
                            .long("disable-illumina")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Disable illumina file naming and read header format. Output file names and reads' header using MGI format."
                            )
                    )
                    .arg(
                        Arg::new("arg_writing_buffer_size")
                            .long("writing-buffer-size")
                            .default_value("67108864")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The size of the buffer for each sample to be filled with data then written once to the disk."
                            )
                    )
                    .arg(
                        Arg::new("arg_lane")
                            .long("lane")
                            .default_value("")
                            .help("The lane number, required for Illumina format.")
                    )
                    .arg(
                        Arg::new("arg_instrument")
                            .long("instrument")
                            .default_value("")
                            .help("The Id of the instrument required for Illumina format.")
                    )
                    .arg(
                        Arg::new("arg_run")
                            .long("run")
                            .default_value("")
                            .help("The run number, required for Illumina format.")
                    )
                    .arg(
                        Arg::new("arg_force")
                            .long("force")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help(
                                "Force running the tool and overwrite existed output/report files."
                            )
                    )
                    .arg(
                        Arg::new("arg_info_file")
                            .long("info-file")
                            .default_value("")
                            .help(
                                "The path to the info file that contains the run information (similar to `BioInfo.csv` generated by MGI machines under the lane directory). Check the documentation for more details."
                            )
                    )
                    .arg(
                        Arg::new("arg_report_level")
                            .long("report-level")
                            .default_value("2")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The level of reporting. 0 no reports will be generated!, 1 data quality and demultipexing reports. 2: all reports (reports on data quality, demultipexing, undetermined and ambigouse barcodes)."
                            )
                    )
                    .arg(
                        Arg::new("arg_compression_level")
                            .long("compression-level")
                            .default_value("1")
                            .value_parser(clap::value_parser!(u32))
                            .help(
                                "The level of compression (between 0 and 12). 0 is fast but no compression, 12 is slow but high compression."
                            )
                    )
                    .arg(
                        Arg::new("arg_compression_buffer_size")
                            .long("compression-buffer-size")
                            .default_value("131072")
                            .value_parser(clap::value_parser!(usize))
                            .help("The size of the buffer for data compression for each sample.")
                    )
                    .arg(
                        Arg::new("arg_memory")
                            .long("memory")
                            .default_value("0")
                            .value_parser(clap::value_parser!(f64))
                            .help(
                                "The requested maximum memory to be used (in giga byte). Check the documentation for memory optimisation options. Default is 0 then the tool will use the available memory on the machine."
                            )
                    )
                    .arg(
                        Arg::new("arg_check_content")
                            .long("validate")
                            .action(ArgAction::SetTrue)
                            .default_value("false")
                            .help("Validate the content of the fastq files.")
                    )
                    .arg(
                        Arg::new("arg_umi_length")
                            .long("umi-length")
                            .default_value("0")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The length of UMI expected at the end of the read (r1 for single-end, or r2 for paired-end)."
                            )
                    )
                    .arg(
                        Arg::new("arg_sample_index")
                            .long("sample-index")
                            .default_value("1")
                            .value_parser(clap::value_parser!(usize))
                            .help(
                                "The index of the sample in the sample sheet, needed for file naming."
                            )
                    )
                    .arg(
                        Arg::new("arg_barcode")
                            .long("barcode")
                            .default_value("")
                            .help(
                                "The barcode of the specific sample to calulate the mismatches for the reports. If not provided, no mismtahces will be calculated."
                            )
                    )
                    .arg(
                        Arg::new("arg_sample_label")
                            .long("sample-label")
                            .default_value("")
                            .help(
                                "Sample label to be used for file naming. if not provided, it will be taken from the name of the input file"
                            )
                    )
            )
            .arg_required_else_help(true)
            .get_matches();

        let start = Instant::now();
        let mut log_level = matches
            .get_one::<String>("arg_log_level")
            .unwrap()
            .to_string()
            .to_lowercase();

        if !["error", "warn", "info", "debug", "trace"].contains(&&log_level[..]) {
            log_level = String::from("info");
            warn!(
                "Log level must be from  [error, warn, info, debug, trace]! Default log level value (info) is used for invalid input!"
            );
        }

        Builder::new()
            .filter_level(LevelFilter::from_str(&log_level).unwrap())
            .target(Target::Stdout)
            .init();
        let args: Vec<String> = env::args().collect();
        info!("Complete Command: {}", args.join(" "));
        info!("Exection start time: {:?}", chrono::offset::Local::now());
        info!("Log level is: {}.", log_level);

        match matches.subcommand() {
            Some(("demultiplex", demultiplex_command)) => {
                initiate_demultiplexing(demultiplex_command);
            }
            Some(("report", report_command)) => {
                let arg_ouput_dir: &String = report_command
                    .get_one::<String>("arg_ouput_dir")
                    .unwrap();
                let arg_qc_report_path: Vec<String> = report_command
                    .get_many::<String>("arg_qc_report_path")
                    .unwrap()
                    .map(|it: &String| it.to_string())
                    .collect::<Vec<_>>();
                let arg_lane: &String = report_command.get_one::<String>("arg_lane").unwrap();
                let arg_prefix: &String = report_command.get_one::<String>("arg_prefix").unwrap();

                merge_qc_reports(&arg_qc_report_path, arg_ouput_dir, &arg_lane, &arg_prefix);
            }
            Some(("template", template_command)) => {
                detect_template(template_command);
            }
            Some(("reformat", reformat_command)) => {
                //let arg_input_folder_path: &String = reformat_command.get_one::<String>("arg_input_folder_path").unwrap();
                match reformat(reformat_command) {
                    Ok(_) => {}
                    Err(err) => eprintln!("Error: {}", err),
                };
            }
            Some((command_nm, _)) => {
                error!("Unknown command `{}`. Please enter a command to perform from (demultiplex, report, template, or reformat)!", command_nm);
            }
            None => {
                panic!(
                    "Please enter a command to perform from (demultiplex, report, template, or reformat)!"
                );
            }
        }
        let dur = start.elapsed();
        info!("{} seconds for performing the task.", dur.as_secs());
        info!("Exection end time: {:?}", chrono::offset::Local::now());
    }
}
