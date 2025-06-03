use std::collections::HashMap;
use log::{ info, error };
use std::time::Instant;
use std::collections::HashSet;
use crate::{ run_manager::ReadInfo, RunManager, SampleManager, variables::* };
use crate::file_utils::{ get_buf_writer, create_output_file };
use std::path::PathBuf;
use std::io::Write;
use std::io::BufWriter;
use std::fs::File;

pub fn write_general_info_report(
    sample_information: &Vec<Vec<String>>,
    sample_statistics: &Vec<Vec<u64>>,
    kept_samples: &Vec<usize>,
    run: &String,
    lane: &String,
    output_file: &PathBuf,
    execluded_samples: &Vec<usize>
) {
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
    let mut out_str = String::from(
        "#sample general info\nSample ID\tM Clusters\tMb Yield ≥ Q30\t% R1 Yield ≥ Q30\t% R2 Yield ≥ Q30\t% R3 Yield ≥ Q30\t% Perfect Index\n"
    );
    let non_sample = vec!["undetermined".to_string(), "ambiguous".to_string()];
    let mut all_q_30;
    let exclude_non_sample = true;
    let mut outfile;
    let mut curr_sample_stat: Vec<f64>;
    let mut tmp_val;
    let mut sample_id;
    let mut lane_statistics: Vec<f64> = vec![0 as f64, 0 as f64, 0 as f64, 0 as f64, 0 as f64];
    let mut sample_statistics_copy = sample_statistics.clone();
    for sample_id_itr in 0..sample_statistics_copy.len() {
        sample_statistics_copy[sample_id_itr].push(sample_id_itr as u64);
    }
    let mut execluded_reads: f64 = 0.0;

    let mut unique_sample_ids = HashSet::new();
    for sample_inf in sample_information {
        unique_sample_ids.insert(sample_inf[SAMPLE_COLUMN].clone());
    }

    if unique_sample_ids.len() == sample_information.len() {
        sample_statistics_copy.sort_by(|a, b|
            b[9].cmp(&a[9]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap()))
        );
    }

    for sample_id_itr in 0..sample_statistics_copy.len() {
        sample_id = sample_statistics_copy[sample_id_itr].last().unwrap().clone() as usize;

        if kept_samples.len() > 0 && !kept_samples.contains(&sample_id) {
            continue;
        }
        if execluded_samples.contains(&sample_id) {
            continue;
        }

        curr_sample_stat = sample_statistics_copy[sample_id_itr]
            .iter()
            .map(|x| x.clone() as f64)
            .collect();

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
        if curr_sample_stat[3] == 0.0 {
            tmp_val = 0.0;
        } else {
            tmp_val = curr_sample_stat[0] / curr_sample_stat[3];
        }
        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));
        out_str.push('\t');

        // % R2 > 30
        if curr_sample_stat[4] == 0.0 {
            tmp_val = 0.0;
        } else {
            tmp_val = curr_sample_stat[1] / curr_sample_stat[4];
        }

        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));

        // % R3 > 30
        out_str.push('\t');

        if curr_sample_stat[5] == 0.0 {
            tmp_val = 0.0;
        } else {
            tmp_val = curr_sample_stat[2] / curr_sample_stat[5];
        }
        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));
        out_str.push('\t');

        // Perfect index
        if curr_sample_stat[9] == 0.0 {
            tmp_val = 0.0;
        } else {
            tmp_val = curr_sample_stat[10] / curr_sample_stat[9];
        }

        out_str.push_str(&format!("{:.3?}", tmp_val * 100.0));

        lane_statistics[3] += curr_sample_stat[6] + curr_sample_stat[7];

        curr_sample_id = curr_sample_id.to_lowercase();
        if non_sample.contains(&curr_sample_id) {
            execluded_reads += curr_sample_stat[9] as f64;
        }

        if !exclude_non_sample || !non_sample.contains(&curr_sample_id) {
            lane_statistics[4] += curr_sample_stat[10];
            //println!("{}  ->  {}   {}    {}    {}", sample_information[sample_id][SAMPLE_COLUMN], lane_statistics[4], curr_sample_stat[10], lane_statistics[1], execluded_reads);
        }

        out_str.push('\n');
    }

    //output_file_path = Path::new(&output_file);
    let mut final_out_str = String::from(
        "#Lane statistics\nRun ID-Lane\tMb Total Yield\tM Total Clusters\t% bases ≥ Q30\tMean Quality\t% Perfect Index\n"
    );
    final_out_str.push_str(&run);
    final_out_str.push('-');
    final_out_str.push_str(&lane);
    final_out_str.push('\t');

    final_out_str.push_str(&format!("{}", lane_statistics[0] / 1000000.0));
    final_out_str.push('\t');

    final_out_str.push_str(&format!("{}", lane_statistics[1] / 1000000.0));
    final_out_str.push('\t');

    final_out_str.push_str(&format!("{:.3?}", (lane_statistics[2] / lane_statistics[0]) * 100.0));
    final_out_str.push('\t');
    final_out_str.push_str(&format!("{:.3?}", lane_statistics[3] / lane_statistics[0]));
    final_out_str.push('\t');
    final_out_str.push_str(
        &format!("{:.3?}", (lane_statistics[4] / (lane_statistics[1] - execluded_reads)) * 100.0)
    );
    final_out_str.push('\n');

    final_out_str.push_str(&out_str);

    outfile = create_output_file(output_file);
    outfile.write_all(&final_out_str.as_bytes()).unwrap();
}

pub fn write_index_info_report(
    sample_information: &Vec<Vec<String>>,
    sample_mismatches: &Vec<Vec<u64>>,
    kept_samples: &Vec<usize>,
    max_mismatches: usize,
    output_file: &PathBuf,
    execluded_samples: &Vec<usize>
) {
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

    if unique_sample_ids.len() == sample_information.len() {
        sample_mismatches_copy.sort_by(|a, b|
            b[0].cmp(&a[0]).then_with(|| a.last().unwrap().cmp(&b.last().unwrap()))
        );
    }
    //sample_mismatches_copy.sort_by(|a, b| b[0].cmp(&a[0]));
    for sample_id_itr in 0..sample_mismatches_copy.len() {
        let sample_index = *sample_mismatches_copy[sample_id_itr].last().unwrap() as usize;
        if kept_samples.len() > 0 && !kept_samples.contains(&sample_index) {
            continue;
        }
        if execluded_samples.contains(&sample_id_itr) {
            continue;
        }
        report_str.push_str(&sample_information[sample_index][SAMPLE_COLUMN]);
        for cnt in 1..max_mismatches + 1 {
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
    let mut outfile = create_output_file(&output_file);
    outfile.write_all(&report_str.as_bytes()).unwrap();
}

#[derive(Default)]
pub struct ReportManager {
    total_samples: usize,
    sample_mismatches: Vec<Vec<u64>>,
    sample_statistics: Vec<Vec<u64>>,
    undetermined_barcodes: HashMap<String, u64>,
    ambiguous_barcodes: HashMap<String, u64>,
}

impl ReportManager {
    pub fn new(total_samples: usize, allowed_mismatches: usize) -> Self {
        let mut sample_mismatches: Vec<Vec<u64>> = Vec::new();
        let mut sample_statistics: Vec<Vec<u64>> = Vec::new();
        for _ in 0..total_samples {
            sample_mismatches.push(vec![0; 2 * allowed_mismatches + 2]);
            sample_statistics.push(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
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
        }
        Self {
            total_samples,
            sample_mismatches,
            sample_statistics,
            undetermined_barcodes: HashMap::new(),
            ambiguous_barcodes: HashMap::new(),
        }
    }

    pub fn update_stats(&mut self, sample_id: usize, index: usize, increment: u64) {
        /*debug!(
            "update stats: sample_id: {}, stat_id: {}, increment: {}",
            sample_id,
            index,
            increment
        );*/
        self.sample_statistics[sample_id][index] += increment;
    }

    pub fn add_stats_entry(&mut self) {
        self.sample_statistics.push(vec![0; 10]);
    }
    pub fn add_mismatches_entry(&mut self, size: usize) {
        self.sample_mismatches.push(vec![0; size]);
    }
    pub fn truncate_mismatches(&mut self, size: usize) {
        for i in 0..self.sample_mismatches.len() {
            self.sample_mismatches[i].truncate(size);
        }
    }

    pub fn update_mismatches(&mut self, sample_id: usize, index: usize, increment: u64) {
        /*debug!(
            "update mismatches: sample_id: {}, mismatches: {}, increment: {}",
            sample_id,
            index,
            increment
        );
        debug!("update mismatches: {:?}", &self.sample_mismatches[sample_id]);
        */ self.sample_mismatches[sample_id][index] += increment;
    }

    pub fn update_undetermined(&mut self, curr_barcode: String, increment: u64) {
        match self.undetermined_barcodes.get_mut(&curr_barcode) {
            Some(barcode_info) => {
                *barcode_info += increment;
            }
            None => {
                self.undetermined_barcodes.insert(curr_barcode, increment);
            }
        };
    }

    pub fn update_ambiguous(&mut self, curr_barcode: String, increment: u64) {
        match self.ambiguous_barcodes.get_mut(&curr_barcode) {
            Some(barcode_info) => {
                *barcode_info += increment;
            }
            None => {
                self.ambiguous_barcodes.insert(curr_barcode, increment);
            }
        };
    }

    pub fn get_sample_reads(&self, sample_id: usize) -> u64 {
        self.sample_mismatches[sample_id][0]
    }

    pub fn get_total_reads(&self) -> u64 {
        let mut reads: u64 = 0;
        for i in 0..self.total_samples {
            for j in 1..self.sample_mismatches[i].len() {
                reads += self.sample_mismatches[i][j];
            }
        }
        reads
    }

    pub fn set_sample_total_reads(&mut self) {
        for i in 0..self.total_samples {
            let mut reads: u64 = 0;
            for j in 1..self.sample_mismatches[i].len() {
                reads += self.sample_mismatches[i][j];
            }
            self.sample_mismatches[i][0] = reads;
        }
    }

    pub fn update(&mut self, report_manager: &ReportManager) {
        if self.sample_mismatches.len() != report_manager.sample_mismatches.len() {
            error!(
                "update ReportManager expects source and destination to have same dimentions. found {} vs {}!",
                self.sample_mismatches.len(),
                report_manager.sample_mismatches.len()
            );
        }
        for i in 0..self.sample_mismatches.len() {
            if self.sample_mismatches[i].len() != report_manager.sample_mismatches[i].len() {
                error!(
                    "update ReportManager expects source and destination to have same dimentions. found {} vs {}!",
                    self.sample_mismatches.len(),
                    report_manager.sample_mismatches.len()
                );
            }
            for j in 0..self.sample_mismatches[i].len() {
                self.sample_mismatches[i][j] += report_manager.sample_mismatches[i][j];
            }
        }

        if self.sample_statistics.len() != report_manager.sample_statistics.len() {
            error!(
                "update ReportManager expects source and destination to have same dimentions. found {} vs {}!",
                self.sample_statistics.len(),
                report_manager.sample_statistics.len()
            );
        }
        for i in 0..self.sample_statistics.len() {
            if self.sample_statistics[i].len() != report_manager.sample_statistics[i].len() {
                error!(
                    "update ReportManager expects source and destination to have same dimentions. found {} vs {}!",
                    self.sample_statistics.len(),
                    report_manager.sample_statistics.len()
                );
            }
            for j in 0..self.sample_statistics[i].len() {
                self.sample_statistics[i][j] += report_manager.sample_statistics[i][j];
            }
        }

        for (key, value) in &report_manager.ambiguous_barcodes {
            self.ambiguous_barcodes
                .entry(key.clone())
                .and_modify(|v| {
                    *v += value;
                })
                .or_insert(*value);
        }

        for (key, value) in &report_manager.undetermined_barcodes {
            self.undetermined_barcodes
                .entry(key.clone())
                .and_modify(|v| {
                    *v += value;
                })
                .or_insert(*value);
        }
    }

    pub fn prepare_final_data(
        &mut self,
        max_mismatches: usize,
        shift: usize,
        barcode_read_info: &ReadInfo,
        paired_read_info: &ReadInfo,
        barcode_length: usize
    ) {
        //let max_mismatches = if all_index_error {allowed_mismatches + 1} else {allowed_mismatches * 2 + 1};
        for sample_id in 0..self.total_samples {
            self.sample_statistics[sample_id][3] =
                (*paired_read_info.sequence_length() as u64) * self.sample_mismatches[sample_id][0];
            self.sample_statistics[sample_id][6] -= self.sample_statistics[sample_id][3] * 33;
            self.sample_statistics[sample_id][shift + 3] =
                ((barcode_read_info.sequence_length() - barcode_length) as u64) *
                self.sample_mismatches[sample_id][0];
            self.sample_statistics[sample_id][5] =
                (barcode_length as u64) * self.sample_mismatches[sample_id][0];
            self.sample_statistics[sample_id][6 + shift] -=
                self.sample_statistics[sample_id][shift + 3] * 33;
            self.sample_statistics[sample_id][8] -= self.sample_statistics[sample_id][5] * 33;
            self.sample_statistics[sample_id][9] = self.sample_mismatches[sample_id][0];
            for cnt in 1..max_mismatches + 1 {
                self.sample_statistics[sample_id].push(self.sample_mismatches[sample_id][cnt]);
            }
        }
    }

    pub fn write_reports(
        &self,
        run_manager: &RunManager,
        sample_manager: &SampleManager,
        reporting_level: usize,
        report_limit: usize,
        mut max_mismatches: usize,
        individual_sample: usize
    ) {
        let start_logs = Instant::now();
        let mut sample_stats_width: usize = self.sample_statistics[0].len();
        let mut execluded_samples = if individual_sample > self.total_samples {
            Vec::new()
        } else {
            let mut execluded_samples: Vec<usize> = (0..self.total_samples).collect();
            execluded_samples.remove(individual_sample);
            execluded_samples
        };

        //println!("{:?}", execluded_samples);
        //println!("Mism: {:?}", self.sample_mismatches);
        //println!("Stat: {:?}", self.sample_statistics);

        if self.sample_mismatches[self.total_samples - 2][0] == 0 {
            execluded_samples.push(self.total_samples - 2);
        }
        if self.sample_mismatches[self.total_samples - 1][0] == 0 {
            execluded_samples.push(self.total_samples - 1);
        }

        let mut report_path_main = if individual_sample == usize::MAX {
            let mut report_path_main = String::from(run_manager.flowcell());
            report_path_main.push('.');
            report_path_main.push_str(run_manager.lane());
            report_path_main
        } else {
            for i in (11..self.sample_statistics[0].len()).rev() {
                if self.sample_statistics[0][i] == 0 {
                    max_mismatches = i - 10;
                    sample_stats_width = i;
                } else {
                    break;
                }
            }
            String::from(sample_manager.sample_information()[0][SAMPLE_COLUMN].clone())
        };
        report_path_main.push_str(&".mgikit.");
        let mut outfile;

        let mut file_name_extra: String;
        let project_samples = sample_manager.project_samples();
        let sample_information = sample_manager.sample_information();
        if individual_sample == usize::MAX {
            for (project_id, samples) in project_samples.iter() {
                if project_id != "." {
                    info!(
                        "generating report for job_number: {} with {} samples.",
                        project_id,
                        samples.len()
                    );
                    file_name_extra = project_id.clone();
                    file_name_extra.push('_');
                } else {
                    info!(
                        "generating report for the whole run with {} samples.",
                        self.sample_mismatches.len()
                    );
                    file_name_extra = String::new();
                }
                let out_file = &run_manager
                    .report_dir()
                    .clone()
                    .join(format!("{}{}info", &file_name_extra, &report_path_main));
                write_index_info_report(
                    sample_information,
                    &self.sample_mismatches,
                    samples,
                    max_mismatches,
                    &out_file,
                    &execluded_samples
                );
                //Finish writing info report

                if reporting_level > 0 {
                    let out_file = &run_manager
                        .report_dir()
                        .clone()
                        .join(format!("{}{}general", &file_name_extra, &report_path_main));

                    write_general_info_report(
                        sample_information,
                        &self.sample_statistics,
                        samples,
                        &run_manager.flowcell(),
                        &run_manager.lane(),
                        &out_file,
                        &execluded_samples
                    );
                }
                //start writing general report
            }
        }
        if reporting_level > 0 {
            let mut out_str = String::from(
                format!(
                    "job_number\tsample_id\tr1_qc_30\tr2_qc_30\tr3_qc_30\tr1_bases\tr2_bases\tr3_bases\tr1_qc\tr2_qc\tr3_qc\tall_reads"
                )
            );
            for cnt in 0..max_mismatches {
                out_str.push('\t');
                out_str.push_str(&cnt.to_string());
                out_str.push_str(&"-mismatches");
            }
            out_str.push('\n');
            for sample_index in 0..self.sample_statistics.len() {
                if !execluded_samples.contains(&sample_index) {
                    out_str.push_str(&sample_information[sample_index][PROJECT_ID_COLUMN]);
                    out_str.push('\t');
                    out_str.push_str(&sample_information[sample_index][SAMPLE_COLUMN]);
                    for cnt in 0..sample_stats_width {
                        out_str.push('\t');
                        out_str.push_str(&self.sample_statistics[sample_index][cnt].to_string());
                    }
                    out_str.push('\n');
                }
            }
            outfile = create_output_file(
                &run_manager.report_dir().clone().join(format!("{}sample_stats", &report_path_main))
            );
            outfile.write_all(&out_str.as_bytes()).unwrap();
        }

        if reporting_level > 1 {
            let mut outfile: BufWriter<File>;
            let mut rep_itr = 0;
            let mut ambiguous_barcodes_out: Vec<_> = self.ambiguous_barcodes.iter().collect();
            if ambiguous_barcodes_out.len() > 0 {
                outfile = get_buf_writer(
                    &run_manager
                        .report_dir()
                        .clone()
                        .join(format!("{}ambiguous_barcode", &report_path_main))
                );
                ambiguous_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
                rep_itr = 0;
                for barcode in &ambiguous_barcodes_out {
                    outfile
                        .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                        .unwrap();
                    rep_itr += 1;
                    if rep_itr == report_limit {
                        break;
                    }
                }

                outfile = get_buf_writer(
                    &run_manager
                        .report_dir()
                        .clone()
                        .join(format!("{}ambiguous_barcode.complete", &report_path_main))
                );
                for barcode in &ambiguous_barcodes_out {
                    outfile
                        .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                        .unwrap();
                }
            }

            //return;

            let mut undetermined_barcodes_out: Vec<_> = self.undetermined_barcodes.iter().collect();
            if undetermined_barcodes_out.len() > 0 {
                outfile = get_buf_writer(
                    &run_manager
                        .report_dir()
                        .clone()
                        .join(&format!("{}undetermined_barcode", report_path_main))
                );
                undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
                for barcode in &undetermined_barcodes_out {
                    outfile
                        .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                        .unwrap();
                    rep_itr += 1;
                    if rep_itr == report_limit {
                        break;
                    }
                }

                outfile = get_buf_writer(
                    &run_manager
                        .report_dir()
                        .clone()
                        .join(&format!("{}undetermined_barcode.complete", report_path_main))
                );
                undetermined_barcodes_out.sort_by(|a, b| (a.1, a.0).cmp(&(b.1, b.0)).reverse());
                for barcode in &undetermined_barcodes_out {
                    outfile
                        .write_all(&format!("{}\t{}\n", barcode.0, barcode.1).as_bytes())
                        .unwrap();
                }
            }
        }

        let log_dur = start_logs.elapsed();

        info!("Writing all logs and reports took {} secs.", log_dur.as_secs());
    }
}
