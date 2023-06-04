
use crate::variables::*;
use crate::sequence_utils::*;
use std::io;
use std::fs;
use std::collections::{HashMap, HashSet};
use std::path::Path;

pub fn parse_template(template: &String) -> Result<[usize; 10], &'static str> {
    let mut template_ls = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    //let vals: Vec<&str> = template.rsplit(':').unwrap();
    let mut shift = 0;
    let mut length;
    //println!("template: {}", template);
    for templ_item in template.rsplit(':') {
        length = templ_item[2..].parse::<usize>().unwrap();
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

pub fn parse_sample_index(
    filename: &Path,
    template: &String,
    i7_rc: bool,
    i5_rc: bool,
) -> Result<
    (
    Vec<(
        u32,
        HashSet<String>,
        HashSet<String>,
        String,
        HashMap<String, (usize, HashMap<String, usize>)>,
        bool,
        [usize; 10],
    )>, 
    Vec<Vec<String>>,
    HashMap<String, Vec<usize>>
    ),
    io::Error,
> {
    //println!("{} -> {}", template, filename.display());

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
        ),
    > = HashMap::new();

    let mut sample_information: Vec<Vec<String>> = Vec::new();
    let mut header: Vec<String> = Vec::new();
    let file_content = fs::read_to_string(filename)?;
    let lines = file_content.lines();
    let mut sample_itr = 0;
    let mut curr_sample_id = usize::MAX;
    let mut curr_template_id = usize::MAX;
    let mut curr_i7 = usize::MAX;
    let mut curr_i5 = usize::MAX;
    let mut curr_i7_rc = usize::MAX;
    let mut curr_i5_rc = usize::MAX;
    let mut curr_project_id = usize::MAX;
    let mut barcode_length: usize = 0;
    let mut curr_sample_info : Vec<String> ; //= Vec::with_capacity(7);
    let mut project_samples : HashMap<String, Vec<usize>> = HashMap::new();
    project_samples.insert(".".to_string(), Vec::new().to_owned());
    let mut report_unassigned_samples = false;
    for line in lines {
        //println!("{}", line);
        if line.len() < 5 {
            continue;
        }

        if header.len() == 0{
            header = line.to_lowercase().split('\t').map(|x| x.trim().to_string()).collect();
            if template == "" && ! header.contains(&String::from("template")) { 
                panic!("Template should be provided either as a general template or within the index/sample file!");
            }
            for header_itr in 0..header.len(){
                if header[header_itr] == "sample_id"{
                    curr_sample_id = header_itr;
                }else if header[header_itr] == "template"{
                    curr_template_id = header_itr;
                }else if header[header_itr] == "i7"{
                    curr_i7 = header_itr;
                }else if header[header_itr] == "i5"{
                    curr_i5 = header_itr;
                }else if header[header_itr] == "i7_rc"{
                    curr_i7_rc = header_itr;
                }else if header[header_itr] == "i5_rc"{
                    curr_i5_rc = header_itr;
                }else if header[header_itr] == "job_number"{
                    curr_project_id = header_itr;
                }                  
            }
        }else{
            let vals: Vec<String> = line.split('\t').map(|x| x.trim().to_string()).collect();
            curr_sample_info = Vec::with_capacity(7) ;

            if curr_sample_id == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_sample_id].to_string());
            }
            
            if curr_i7 == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_i7].to_string());
            }

            if curr_i5 == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_i5].to_string());
            }
            
            if curr_template_id == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_template_id].to_string());
            }

            if curr_i7_rc == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_i7_rc].to_string());
            }

            if curr_i5_rc == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_i5_rc].to_string());
            }
            
            if curr_project_id == usize::MAX {
                curr_sample_info.push(".".to_string());
            }else{
                curr_sample_info.push(vals[curr_project_id].to_string());
            }
            
            match project_samples.get_mut(&curr_sample_info[PROJECT_ID_COLUMN]) {
                Some(curr_project_samples) => {curr_project_samples.push(sample_itr);},
                None => {
                    project_samples.insert(curr_sample_info[PROJECT_ID_COLUMN].to_owned(), vec!(sample_itr).to_owned());
                }
            };
            
            /*
            println!("{} ***  {:?}   *** {}",curr_sample_info[PROJECT_ID_COLUMN], project_samples["."], curr_sample_info[PROJECT_ID_COLUMN] != "." && project_samples["."].len() > 0 );
            if curr_sample_info[PROJECT_ID_COLUMN] != "." && project_samples["."].len() > 0{
                panic!("All samples should have a project_id or none of them!");
            }
            */

            if ! report_unassigned_samples && project_samples.len() > 1 && project_samples["."].len() > 0{
                println!("Some samples are not assigned to any project, they will be in the reports of the whole run!");
                report_unassigned_samples = true;
            }
        
            let curr_template = if template == "" {
                curr_sample_info[TEMPLATE_COLUMN].clone()
            } else {
                template.clone()
            };
            //println!("{:?}", curr_sample_info);
            //println!("curr {}   tmpl: {}",curr_template_id, curr_template);

            if curr_sample_info[I7_COLUMN] == "." || curr_sample_info[I7_COLUMN].len() < 3 {
                panic!("i7 ({}) should be longer than 3 chars!", vals[I7_COLUMN]);
            }

            if (curr_sample_info[I5_COLUMN] == "." || curr_sample_info[I5_COLUMN].len() < 3) && curr_template.contains("i5") {
                panic!(
                    "i5 ({}) should be longer than 3 chars! or the template should not contains i5",
                    curr_sample_info[I5_COLUMN]
                );
            }

            let mut check_i5 = false;
            if curr_sample_info[I5_COLUMN] != "." {
                check_i5 = true;
            }

            let i7 = if (template.len() > 0 && i7_rc)
                || (template.len() == 0 && curr_sample_info[I7_RC_COLUMN] == "1")
            {
                //println!("revi7 {}  {}", template, i7_rc);
                reverse_complement(&curr_sample_info[I7_COLUMN]).unwrap()
            } else {
                curr_sample_info[I7_COLUMN].clone()
            };
            //println!("I7s: {} -> {}", vals[I7_COLUMN], i7);

            let i5 = if check_i5
                && ((template.len() > 0 && i5_rc) || (template.len() == 0 && curr_sample_info[I5_RC_COLUMN] == "1"))
            {
                reverse_complement(&curr_sample_info[I5_COLUMN]).unwrap()
            } else if check_i5 {
                curr_sample_info[I5_COLUMN].clone()
            } else {
                String::new()
            };

            //println!("II55  {}", i5);
            match template_ls.get_mut(&curr_template) {
                Some(tmp_data) => {
                    tmp_data.0 += 1;
                    tmp_data.1.insert(i7.clone());
                    tmp_data.2.insert(i5.clone());

                    match tmp_data.4.get_mut(&i7) {
                        Some(i7_item) => {
                            if check_i5 {
                                match i7_item.1.get_mut(&i5) {
                                    Some(_) => panic!("Two samples having the same indexes! i7: {} and i5: {}", &i7, &i5),
                                    None => {
                                        i7_item.1.insert(i5, sample_itr);
                                    }
                                }
                            } else {
                                panic!("Two samples having the same i7 indexes! i7: {}", &i7)
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
                        if tmp.6[9] != barcode_length{
                            panic!("The barcode length should be the same for all samples!");
                        }
                    }else{
                        barcode_length = tmp.6[9];
                    }
                    template_ls.insert(curr_template, tmp);
                }
            };

            sample_information.push(curr_sample_info.to_owned());
            //break;
            //println!("After_itr {:?}", template_ls.get(&curr_tmpllll).unwrap().4.get(&i77).unwrap().1.get(&i55).unwrap());

            //output.push(vals.to_owned());
            //println!("line: {}", vals.join(";"));
            
            sample_itr += 1;
        }
        
    }

    if sample_information.len() == 0{
        panic!("Sample sheet seems to be empty! No sample is found!");
    }
    //println!("final {:?}", template_ls.get(&String::from("i78:i58")).unwrap().4.get(&String::from("CATGCCTA")).unwrap().1.get(&String::from("ATAGAGAG")).unwrap());

    let mut out_template_data: Vec<_> = template_ls.values().map(|x| x.to_owned()).collect();
    if out_template_data.len() > 1 {
        out_template_data.sort_by(|a, b| a.0.cmp(&b.0).reverse());
    }

    /*
    for v in &out_template_data {
    println!("tmpl {} ->  cnount: {}\ni7s: {:?}\ni5s: {:?}\nchecki5 {} ", v.3, v.0, v.1, v.2, v.5);

    for (k7, v7) in &v.4{
    if !v.5{
    println!("only I7: {} {:?}", k7, v7.0);
    }else{
    for (k75, v5) in &v7.1{
    println!("{} and {} -> {:?}", k7, k75, v5);
    }
    }
    println!("\n----------------------\n");
    }
    println!("\n******************************\n");

    }
    */

    //println! ("{:?}\n*********\n{:?}", out_template_data, sample_information);
    //panic!("done!");
    match project_samples.get_mut(&String::from(".")) {
        Some(curr_project_samples) => {curr_project_samples.clear();},
        None => {
            panic!("There must be a project for the whole run with the '.' name!");
        }
    };

    Ok((out_template_data, 
        sample_information,
        project_samples))
}

