use getset::Getters;
use crate::variables::*;
use crate::file_utils::*;
use core::panic;
use std::io;
use std::fs;
use std::collections::{ HashMap, HashSet };
use std::path::Path;
use std::path::PathBuf;
use log::info;
use itertools::Itertools;

#[derive(Getters, Clone, Default)]
pub struct SampleManager {
    #[getset(get = "pub")]
    sample_information: Vec<Vec<String>>,
    #[getset(get = "pub")]
    project_samples: HashMap<String, Vec<usize>>,
    #[getset(get = "pub")]
    writing_samples: Vec<usize>,
    #[getset(get = "pub")]
    unique_samples_ids: Vec<usize>,
    #[getset(get = "pub")]
    all_template_data: Vec<
        (
            u32,
            HashSet<String>,
            HashSet<String>,
            String,
            HashMap<String, (usize, HashMap<String, usize>)>,
            bool,
            [usize; 10],
        )
    >,
}

impl SampleManager {
    pub fn new<P: Into<PathBuf>>(
        sample_sheet_path: P,
        template: String,
        i7_rc: bool,
        i5_rc: bool,
        undetermined_label: String,
        ambiguous_label: String
    ) -> Self {
        let sample_sheet_path = sample_sheet_path.into();
        if !sample_sheet_path.exists() {
            panic!("Sample sheet file is invalid!");
        }
        check_file(&sample_sheet_path);

        let mut sample_information = load_sample_sheet(&sample_sheet_path).unwrap();
        info!("{} Samples were found in the input sample sheet.", sample_information.len());

        let (w_s, u_s) = get_writing_unique_samples(&sample_information).unwrap();
        // parse sample/index file and get all mismatches
        if template.len() > 0 {
            info!("General template is provided and will be used for all samples: {}", &template);
            if i7_rc {
                info!("i7 will be converted to the reverse complement!");
            }

            if i5_rc {
                info!("i5 will be converted to the reverse complement!");
            }
        } else {
            info!("Template will be used from sample/index map file.");
            //i7_rc = false;
            //i5_rc = false;
        }

        let all_template_data = extract_templates_information(
            &sample_information,
            &template,
            i7_rc,
            i5_rc
        ).unwrap();
        if all_template_data.len() > 1 {
            info!("Mixed library is detected! different barcode templates for some samples!");
        } else {
            info!("Same barcode template is used for all samples!");
        }
        sample_information.push(
            vec![
                undetermined_label.clone(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::from(".")
            ]
        );
        sample_information.push(
            vec![
                ambiguous_label.clone(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::from(".")
            ]
        );

        Self {
            project_samples: extract_project_samples(&sample_information).unwrap(),
            all_template_data,
            writing_samples: w_s,
            unique_samples_ids: u_s,
            sample_information,
        }
    }

    pub fn from_simple_sheet<P: Into<PathBuf>>(sample_sheet_path: P) -> Self {
        let sample_sheet_path = sample_sheet_path.into();
        if !sample_sheet_path.exists() {
            panic!("Sample sheet file is invalid!");
        }
        check_file(&sample_sheet_path);

        let sample_information = load_sample_sheet(&sample_sheet_path).unwrap();
        info!("{} Samples were found in the input sample sheet.", sample_information.len());
        Self {
            project_samples: extract_project_samples(&sample_information).unwrap(),
            all_template_data: Vec::new(),
            writing_samples: Vec::new(),
            unique_samples_ids: Vec::new(),
            sample_information,
        }
    }

    pub fn dummy_sample(sample_label: String) -> Self {
        let mut sample_information = Vec::new();
        sample_information.push(
            vec![
                sample_label.clone(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::from(".")
            ]
        );
        sample_information.push(
            vec![
                String::from("Undetermined"),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::from(".")
            ]
        );
        sample_information.push(
            vec![
                String::from("Ambiguous"),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::from(".")
            ]
        );

        let (w_s, u_s) = get_writing_unique_samples(&sample_information).unwrap();
        // parse sample/index file and get all mismatches
        Self {
            project_samples: extract_project_samples(&sample_information).unwrap(),
            all_template_data: Vec::new(),
            writing_samples: w_s,
            unique_samples_ids: u_s,
            sample_information,
        }
    }

    pub fn add_sample(&mut self, sample_label: &String, project_label: &String) {
        self.sample_information.push(
            vec![
                sample_label.clone(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                project_label.clone()
            ]
        );
    }

    pub fn pop_sample(&mut self) {
        self.sample_information.pop();
    }

    pub fn add_project_sample(&mut self, project: &String, sample_indx: usize) {
        match self.project_samples.get_mut(project) {
            Some(samples) => {
                if !samples.contains(&sample_indx) {
                    samples.push(sample_indx);
                }
            }
            None => {
                self.project_samples.insert(project.clone(), vec![sample_indx]);
            }
        };
    }

    pub fn get_sample_count(&self) -> usize {
        self.sample_information.len().clone()
    }

    pub fn get_sample_index(&self, sample_id: &String) -> usize {
        for i in 0..self.sample_information.len() {
            if self.sample_information[i][0] == *sample_id {
                return i;
            }
        }
        return usize::MAX;
    }

    pub fn get_samples_indices(&self) -> Result<Vec<[String; 4]>, io::Error> {
        let mut sample_indexes: Vec<[String; 4]> = Vec::new();
        for sample_info in self.sample_information() {
            //println!("{:?}", &sample_info);
            if sample_info[I5_COLUMN] != ".".to_string() {
                sample_indexes.push([
                    sample_info[I7_COLUMN].to_string(),
                    reverse_complement(&sample_info[I7_COLUMN]).unwrap(),
                    sample_info[I5_COLUMN].to_string(),
                    reverse_complement(&sample_info[I5_COLUMN]).unwrap(),
                ]);
            } else {
                sample_indexes.push([
                    sample_info[I7_COLUMN].to_string(),
                    reverse_complement(&sample_info[I7_COLUMN]).unwrap(),
                    String::new(),
                    String::new(),
                ]);
            }
        }
        //println!("{:?}", &sample_indexes);
        Ok(sample_indexes)
    }
}

pub fn parse_template(template: &String) -> Result<[usize; 10], &'static str> {
    let mut template_ls = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    //let vals: Vec<&str> = template.rsplit(':').unwrap();
    let mut shift = 0;
    let mut length;
    //debug!("template: {}", template);
    for templ_item in template.rsplit(':') {
        if !["i7", "i5", "um", "--"].contains(&&templ_item[0..2]) {
            panic!("template ({}) does not match the expected format as explianed in mgikit documenation!", template);
        }
        length = templ_item[2..]
            .parse::<usize>()
            .expect(
                &format!("template ({}) does not match the expected format as explianed in mgikit documenation!", template)
            );
        shift = shift + length;
        if templ_item[0..2] == *"i7" {
            template_ls[0] = 1;
            template_ls[1] = length;
            template_ls[2] = shift;
        } else if templ_item[0..2] == *"i5" {
            template_ls[3] = 1;
            template_ls[4] = length;
            template_ls[5] = shift;
        } else if templ_item[0..2] == *"um" {
            template_ls[6] = 1;
            template_ls[7] = length;
            template_ls[8] = shift;
        } else if templ_item[0..2] != *"--" {
            //println!("Should not be here!");
            return Err("Wrong format for template!");
        }
    }
    template_ls[9] = shift;
    /*println!("{} -> {}.{}.{}    {}.{}.{}   {}.{}.{}", template, template_ls[0][0],template_ls[0][1],template_ls[0][2],
    template_ls[1][0],template_ls[1][1],template_ls[1][2],
    template_ls[2][0],template_ls[2][1],template_ls[2][2]);*/

    Ok(template_ls)
}

pub fn load_sample_sheet(filename: &Path) -> Result<Vec<Vec<String>>, io::Error> {
    let mut sample_information: Vec<Vec<String>> = Vec::new();
    let mut header: Vec<String> = Vec::new();
    let file_content = fs::read_to_string(filename)?;
    let lines = file_content.lines();
    let mut curr_sample_id = usize::MAX;
    let mut curr_template_id = usize::MAX;
    let mut curr_i7 = usize::MAX;
    let mut curr_i5 = usize::MAX;
    let mut curr_i7_rc = usize::MAX;
    let mut curr_i5_rc = usize::MAX;
    let mut curr_project_id = usize::MAX;
    let mut curr_sample_info: Vec<String>; //= Vec::with_capacity(7);
    let mut delimiter = '\t';
    for line in lines {
        //println!("ZZZ: {}", line);
        if line.trim().len() < 5 {
            continue;
        }
        if header.len() == 0 {
            header = line
                .to_lowercase()
                .split('\t')
                .map(|x| x.trim().to_string())
                .collect();
            if header.len() < 2 {
                header = line
                    .to_lowercase()
                    .split(',')
                    .map(|x| x.trim().to_string())
                    .collect();
                delimiter = ',';
                if header.len() < 2 {
                    panic!("Sample sheet columns should be separated by ',' or '\t'!");
                }
            }
            for header_itr in 0..header.len() {
                if header[header_itr] == "sample_id" {
                    curr_sample_id = header_itr;
                } else if header[header_itr] == "template" {
                    curr_template_id = header_itr;
                } else if header[header_itr] == "i7" {
                    curr_i7 = header_itr;
                } else if header[header_itr] == "i5" {
                    curr_i5 = header_itr;
                } else if header[header_itr] == "i7_rc" {
                    curr_i7_rc = header_itr;
                } else if header[header_itr] == "i5_rc" {
                    curr_i5_rc = header_itr;
                } else if header[header_itr] == "job_number" {
                    curr_project_id = header_itr;
                }
            }
        } else {
            let vals: Vec<String> = line
                .split(delimiter)
                .map(|x| x.trim().to_string())
                .collect();
            //println!("values of sample: {:?}", vals);
            curr_sample_info = Vec::with_capacity(7);

            if curr_sample_id == usize::MAX || vals[curr_sample_id].to_string().len() == 0 {
                panic!(
                    "sample_id column is mandatory in the samplesheet and must not be an empty string!"
                );
            } else {
                curr_sample_info.push(vals[curr_sample_id].to_string());
            }

            if curr_i7 == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                if
                    !vals[curr_i7]
                        .to_string()
                        .chars()
                        .all(|c| ['A', 'C', 'G', 'T'].contains(&c))
                {
                    panic!(
                        "Index must only contain A, C, G and T letters! found '{}'",
                        vals[curr_i7].to_string()
                    );
                }
                curr_sample_info.push(vals[curr_i7].to_string());
            }

            if curr_i5 == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                curr_sample_info.push(vals[curr_i5].to_string());
            }

            if curr_template_id == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                curr_sample_info.push(vals[curr_template_id].to_string());
            }

            if curr_i7_rc == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                if
                    vals[curr_i7_rc].to_string() != "." &&
                    vals[curr_i7_rc].to_string() != "0" &&
                    vals[curr_i7_rc].to_string() != "1"
                {
                    panic!(
                        "i7_rc must be either '.', '0' or '1' when reverse complementary! found '{}'",
                        vals[curr_i7_rc].to_string()
                    );
                }
                curr_sample_info.push(vals[curr_i7_rc].to_string());
            }

            if curr_i5_rc == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                if
                    vals[curr_i5_rc].to_string() != "." &&
                    vals[curr_i5_rc].to_string() != "0" &&
                    vals[curr_i5_rc].to_string() != "1"
                {
                    panic!(
                        "i5_rc must be either '.', '0' or '1' when reverse complementary! found '{}'",
                        vals[curr_i5_rc].to_string()
                    );
                }
                curr_sample_info.push(vals[curr_i5_rc].to_string());
            }

            if curr_project_id == usize::MAX {
                curr_sample_info.push(".".to_string());
            } else {
                curr_sample_info.push(vals[curr_project_id].to_string());
            }

            if curr_sample_info[I7_COLUMN] == "." || curr_sample_info[I7_COLUMN].len() < 3 {
                panic!("i7 ({}) should be longer than 3 chars!", vals[I7_COLUMN]);
            }

            if curr_sample_info[I5_COLUMN] != "." {
                if curr_sample_info[I5_COLUMN].len() < 3 {
                    panic!("i5 ({}) should be longer than 3 chars!", vals[I5_COLUMN]);
                }
                if
                    !vals[curr_i5]
                        .to_string()
                        .chars()
                        .all(|c| ['A', 'C', 'G', 'T'].contains(&c))
                {
                    panic!(
                        "Index must only contain A, C, G and T letters! found '{}'",
                        vals[curr_i5].to_string()
                    );
                }
            }

            sample_information.push(curr_sample_info.to_owned());
        }
    }

    if sample_information.len() == 0 {
        panic!("Sample sheet seems to be empty! No sample is found!");
    }

    Ok(sample_information)
}

pub fn extract_project_samples(
    sample_information: &Vec<Vec<String>>
) -> Result<HashMap<String, Vec<usize>>, io::Error> {
    //println!("{} -> {}", template, filename.display());

    let mut sample_itr = 0;
    let mut project_samples: HashMap<String, Vec<usize>> = HashMap::new();
    project_samples.insert(".".to_string(), Vec::new().to_owned());
    let mut report_unassigned_samples = false;
    for curr_sample_info in sample_information {
        match project_samples.get_mut(&curr_sample_info[PROJECT_ID_COLUMN]) {
            Some(curr_project_samples) => {
                curr_project_samples.push(sample_itr);
            }
            None => {
                project_samples.insert(
                    curr_sample_info[PROJECT_ID_COLUMN].to_owned(),
                    vec![sample_itr].to_owned()
                );
            }
        }

        if
            !report_unassigned_samples &&
            project_samples.len() > 1 &&
            project_samples["."].len() > 0
        {
            info!(
                "Some samples are not assigned to any project, they will be in the reports of the whole run!"
            );
            report_unassigned_samples = true;
        }

        sample_itr += 1;
    }

    match project_samples.get_mut(&String::from(".")) {
        Some(curr_project_samples) => {
            curr_project_samples.clear();
        }
        None => {
            panic!("There must be a project for the whole run with the '.' name!");
        }
    }

    Ok(project_samples)
}

pub fn extract_templates_information(
    sample_information: &Vec<Vec<String>>,
    template: &String,
    i7_rc: bool,
    i5_rc: bool
) -> Result<
    Vec<
        (
            u32,
            HashSet<String>,
            HashSet<String>,
            String,
            HashMap<String, (usize, HashMap<String, usize>)>,
            bool,
            [usize; 10],
        )
    >,
    io::Error
> {
    /*
    template details = (
                    1,
                    all i7s indexes,
                    all i5s indexes,
                    template string,
                    (i7, (i5, sample index)) if no i5, it will be i7 and sample id,
                    boolean to check i5,
                    template information as where to find the indexes,
                );
     */
    //println!("{} -> {}", template, filename.display());
    //debug!("{:?}", sample_information);
    let mut template_ls: HashMap<
        String,
        (
            u32,
            HashSet<String>,
            HashSet<String>,
            String,
            HashMap<String, (usize, HashMap<String, usize>)>,
            bool,
            [usize; 10],
        )
    > = HashMap::new();

    let mut sample_itr = 0;
    let mut barcode_length: usize = 0;
    for curr_sample_info in sample_information {
        let curr_template = if template == "" {
            curr_sample_info[TEMPLATE_COLUMN].clone()
        } else {
            template.clone()
        };
        //debug!("curr template {}  -> *{}*", curr_template, template);

        if
            (curr_sample_info[I5_COLUMN] == "." || curr_sample_info[I5_COLUMN].len() < 3) &&
            curr_template.contains("i5")
        {
            panic!(
                "i5 ({}) should be longer than 3 chars! or the template should not contains i5",
                curr_sample_info[I5_COLUMN]
            );
        }

        let check_i5 = if curr_sample_info[I5_COLUMN].len() > 1 { true } else { false };

        let i7 = if
            (template.len() > 0 && i7_rc) ||
            (template.len() == 0 && curr_sample_info[I7_RC_COLUMN] == "1")
        {
            //println!("revi7 {}  {}", template, i7_rc);
            reverse_complement(&curr_sample_info[I7_COLUMN]).unwrap()
        } else {
            curr_sample_info[I7_COLUMN].clone()
        };
        //println!("I7s: {} -> {}", vals[I7_COLUMN], i7);

        let i5 = if
            check_i5 &&
            ((template.len() > 0 && i5_rc) ||
                (template.len() == 0 && curr_sample_info[I5_RC_COLUMN] == "1"))
        {
            reverse_complement(&curr_sample_info[I5_COLUMN]).unwrap()
        } else if check_i5 {
            curr_sample_info[I5_COLUMN].clone()
        } else {
            String::new()
        };

        match template_ls.get_mut(&curr_template) {
            Some(tmp_data) => {
                tmp_data.0 += 1;
                tmp_data.1.insert(i7.clone());
                tmp_data.2.insert(i5.clone());

                match tmp_data.4.get_mut(&i7) {
                    Some(i7_item) => {
                        if check_i5 {
                            match i7_item.1.get_mut(&i5) {
                                Some(_) =>
                                    panic!(
                                        "Two samples having the same indexes! i7: {} and i5: {}",
                                        &i7,
                                        &i5
                                    ),
                                None => {
                                    i7_item.1.insert(i5, sample_itr);
                                }
                            }
                        } else {
                            panic!("Two samples having the same i7 indexes! i7: {}", &i7);
                        }
                    }
                    None => {
                        if check_i5 {
                            let mut tmp = HashMap::new();
                            tmp.insert(i5, sample_itr);
                            tmp_data.4.insert(i7, (1000, tmp));
                        } else {
                            tmp_data.4.insert(i7, (sample_itr, HashMap::new()));
                        }
                    }
                };
            }
            None => {
                let mut i7s = HashSet::new();
                let mut i5s = HashSet::new();
                i7s.insert(i7.clone());
                if check_i5 {
                    i5s.insert(i5.clone());
                }
                let mut sample_info = HashMap::new();
                if check_i5 {
                    let mut tmp = HashMap::new();
                    //let zzz = i5.clone();
                    //let zzzz = i7.clone();
                    tmp.insert(i5, sample_itr);
                    sample_info.insert(i7, (1000, tmp));
                    //println!("SSSS {:?}", sample_info.get(&zzzz).unwrap().1.get(&zzz).unwrap());
                } else {
                    sample_info.insert(i7, (sample_itr, HashMap::new()));
                }

                let tmp = (
                    1,
                    i7s,
                    i5s,
                    curr_template.clone(),
                    sample_info,
                    check_i5,
                    parse_template(&curr_template).unwrap(),
                );
                if template_ls.len() > 0 {
                    if tmp.6[9] != barcode_length {
                        panic!("The barcode length should be the same for all samples!");
                    }
                } else {
                    barcode_length = tmp.6[9];
                }
                template_ls.insert(curr_template, tmp);
            }
        }

        sample_itr += 1;
    }

    let mut out_template_data: Vec<_> = template_ls
        .values()
        .map(|x| x.to_owned())
        .collect();
    if out_template_data.len() > 1 {
        out_template_data.sort_by(|a, b| {
            if a.0 == b.0 {
                b.0.cmp(&a.0)
                //b.3.cmp(&a.3)
            } else {
                b.0.cmp(&a.0)
            }
            //a.0.cmp(&b.0).reverse()
        });
    }

    Ok(out_template_data)
}

pub fn get_writing_unique_samples(
    sample_information: &Vec<Vec<String>>
) -> Result<(Vec<usize>, Vec<usize>), io::Error> {
    //println!("{} -> {}", template, filename.display());

    let mut sample_itr = 0;
    let mut sample_duplicate: HashMap<String, usize> = HashMap::new();
    let mut writing_samples: Vec<usize> = Vec::new();
    let mut unique_sample_id: Vec<usize> = Vec::new();
    let mut dup_ids = 0;
    let mut curr_unique_id;
    for curr_sample_info in sample_information {
        {
            let duplicate_sample_id = match sample_duplicate.get(&curr_sample_info[SAMPLE_COLUMN]) {
                Some(dup_id) => {
                    curr_unique_id = unique_sample_id[*dup_id];
                    *dup_id
                }
                None => {
                    sample_duplicate.insert(curr_sample_info[SAMPLE_COLUMN].clone(), sample_itr);
                    curr_unique_id = dup_ids;
                    dup_ids += 1;
                    sample_itr
                }
            };
            writing_samples.push(duplicate_sample_id);
            unique_sample_id.push(curr_unique_id);
            sample_itr += 1;
        }
    }

    //adding undetrmined and ambiguose samples ids.
    writing_samples.push(writing_samples.len());
    unique_sample_id.push(unique_sample_id.len());
    writing_samples.push(writing_samples.len());
    unique_sample_id.push(unique_sample_id.len());

    Ok((writing_samples, unique_sample_id))
}

pub fn get_all_mismatches<'a>(
    all_indexes: &'a HashSet<String>,
    allowed_mismatches: usize
) -> HashMap<Vec<u8>, (Vec<&'a String>, usize)> {
    if all_indexes.len() == 0 {
        return HashMap::new();
    }

    let nec_ls = ['A', 'C', 'T', 'G', 'N'];
    //let mut index_combinations: HashMap<&String, Vec<HashSet<String>>> = HashMap::new();

    let mut final_output: HashMap<String, Vec<Vec<&String>>> = HashMap::new();

    let mut possible_mismatches: HashSet<String>;

    for index in all_indexes {
        //println!("index: {} -> {:?}", index, "Z");//, final_output.keys());
        //index_combinations.insert(&index, Vec::new());
        let mut tmp: Vec<Vec<&String>> = vec![Vec::new(); allowed_mismatches + 1];
        tmp[0].push(&index);
        match final_output.get_mut(index) {
            Some(x) => {
                *x = tmp;
                //println!("ovvoride index");
                //panic!("duplicated index!")
            }
            None => {
                final_output.insert(index.clone(), tmp);
            }
        }

        for mismatch_itr in 1..allowed_mismatches + 1 {
            //println!("Iter: {}", mismatch_itr);
            possible_mismatches = HashSet::new();
            let combination_ls = (0..index.len()).combinations(mismatch_itr);
            for combination in combination_ls {
                //println!("Z: {:?}", combination);

                let mut all_combinations: Vec<String> = nec_ls
                    .iter()
                    .filter(|&i| *i != index.chars().nth(combination[0]).unwrap())
                    .map(|x| x.to_string())
                    .collect();
                //println!("ALL: {:?}", all_combinations);
                {
                    for nec_ind in 1..combination.len() {
                        let ls_nec: Vec<char> = nec_ls
                            .iter()
                            .filter(|&i| *i != index.chars().nth(combination[nec_ind]).unwrap())
                            .cloned()
                            .collect();
                        //println!("inner {} -> {:?}", nec_ind, ls_nec);

                        let mut tmp: Vec<String> = Vec::new();
                        for x in all_combinations {
                            for y in &ls_nec {
                                let mut z = x.clone();
                                z.push(*y);
                                //println!("merge: {} - {}  -> {}", x, y, z);
                                tmp.push(z);
                            }
                        }
                        all_combinations = tmp;
                    }

                    for combination_inner in all_combinations {
                        let mut index_edited = index.clone();

                        for nec_ind in 0..combination.len() {
                            //index_edited[combination[nec_ind]] = combination_inner.chars().nth(nec_ind).unwrap();
                            index_edited.replace_range(
                                combination[nec_ind]..combination[nec_ind] + 1,
                                &combination_inner[nec_ind..nec_ind + 1]
                            );
                        }
                        //println!("{}   {}", index, index_edited);

                        match final_output.get_mut(&index_edited) {
                            Some(x) => {
                                x[mismatch_itr].push(index);
                            }
                            None => {
                                let mut tmp = Vec::new();
                                for mis_itr in 0..allowed_mismatches + 1 {
                                    if mis_itr == mismatch_itr {
                                        tmp.push(vec![index]);
                                    } else {
                                        tmp.push(Vec::new());
                                    }
                                }
                                final_output.insert(index_edited.clone(), tmp);
                            }
                        }
                        possible_mismatches.insert(index_edited);
                    }
                }
            }

            /*match index_combinations.get_mut(&index){
                Some(x) => {x.push(possible_mismatches);},
                None => panic!("Something wrong")
            };*/
        }
    }

    let mut sample_index_map = HashMap::new();
    for (key, value) in final_output {
        for i in 0..allowed_mismatches + 1 {
            if value[i].len() > 0 {
                //let mut tmp_z = value[i].clone();
                //tmp_z.sort_by(|a, b| a.cmp(b));
                sample_index_map.insert(key.as_bytes().to_vec(), (value[i].to_owned(), i));
                break;
            }
        }
    }

    sample_index_map
}

pub fn reverse_complement(seq: &String) -> Result<String, &'static str> {
    //println!("rc: *{}*", seq);
    let mut out_string = String::new();
    for nec in seq.chars().rev() {
        out_string.push(match nec {
            'n' => 'N',
            'N' => 'N',
            'A' => 'T',
            'a' => 'T',
            'T' => 'A',
            't' => 'A',
            'C' => 'G',
            'c' => 'G',
            'G' => 'C',
            'g' => 'C',
            _ => panic!("Wrong neucltide!"),
        });
    }
    Ok(out_string)
}

#[cfg(test)]
mod tests {
    use super::reverse_complement;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(&String::from("ACGTGCGC")).unwrap(), "GCGCACGT");
        assert_eq!(reverse_complement(&String::from("GGGGGGGG")).unwrap(), "CCCCCCCC");
        assert_eq!(reverse_complement(&String::from("ATATATNN")).unwrap(), "NNATATAT");
        assert_eq!(reverse_complement(&String::from("agcagccc")).unwrap(), "GGGCTGCT");
    }
}
