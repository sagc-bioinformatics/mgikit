use crate::variables::*;
use crate::sequence_utils::*;
use core::panic;
use std::io;
use std::fs;
use std::collections::{HashMap, HashSet};
use std::path::Path;


/// It takes a template string as input and parse it.  
/// It return a list of 10 elements incouding i7_exist, i7 length, shift from the begining of the barcode. 
/// The three items are repeasted to i5 then umi. `--` are ignored. last item is the length of the barcode.
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

/// It returns (template_list, sample_info, project_samples)
/// 

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
    HashMap<String, Vec<usize>>,
    Vec<usize>,
    Vec<usize> 
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
    let mut sample_duplicate : HashMap<String, usize> = HashMap::new();
    let mut writing_samples: Vec<usize> = Vec::new();
    let mut unique_sample_id: Vec<usize> = Vec::new();
    let mut dup_ids = 0;
    let mut curr_unique_id;
    let mut delimiter = '\t';
    for line in lines {
        //println!("{}", line);
        if line.len() < 5 {
            continue;
        }
        if header.len() == 0{
            header = line.to_lowercase().split('\t').map(|x| x.trim().to_string()).collect();
            if header.len() < 2{
                header = line.to_lowercase().split(',').map(|x| x.trim().to_string()).collect();
                delimiter = ',';
                if header.len() < 2{
                    panic!("Sample sheet columns should be separated by ',' or '\t'!");
                }
            }
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
            let vals: Vec<String> = line.split(delimiter).map(|x| x.trim().to_string()).collect();
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

            let check_i5 = if curr_sample_info[I5_COLUMN].len() > 1 {
                true
            }else{
                false
            };

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

            let duplicate_sample_id = match sample_duplicate.get(&curr_sample_info[SAMPLE_COLUMN]) {
                Some(dup_id) => {
                    curr_unique_id = unique_sample_id[*dup_id]; 
                    *dup_id
                },
                None => {
                    sample_duplicate.insert(curr_sample_info[SAMPLE_COLUMN].clone(), sample_itr);
                    curr_unique_id = dup_ids; 
                    dup_ids += 1;
                    sample_itr
                }
            };
            writing_samples.push(duplicate_sample_id);
            unique_sample_id.push(curr_unique_id);
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
            
            //if duplicate_sample_id == sample_itr {
                sample_information.push(curr_sample_info.to_owned());
            //}
            
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
        out_template_data.sort_by(|a, b| {
            if a.0 == b.0 {
                b.0.cmp(&a.0)
                //b.3.cmp(&a.3)
            } else {
                b.0.cmp(&a.0)
            }
            //a.0.cmp(&b.0).reverse()
        }
        );
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

    //adding undetrmined and ambiguose samples ids.
    writing_samples.push(writing_samples.len());
    unique_sample_id.push(unique_sample_id.len());
    writing_samples.push(writing_samples.len());
    unique_sample_id.push(unique_sample_id.len());

    Ok((out_template_data, 
        sample_information,
        project_samples,
        writing_samples,
        unique_sample_id))
}



pub fn read_sample_sheet_into_dic(
    filename: &Path) -> Result<(Vec<usize>, Vec<Vec<String>>, Vec<[String; 4]>), io::Error> {
    
    let mut sample_information: (Vec<usize>, Vec<Vec<String>>) = (vec![usize::MAX, usize::MAX, usize::MAX, usize::MAX], Vec::new());
    let mut sample_indexes: Vec<[String; 4]> = Vec::new();
    let mut header: Vec<String> = Vec::new();
    let file_content = fs::read_to_string(filename)?;
    let lines = file_content.lines();
    let mut curr_sample_id = usize::MAX;
    let mut curr_i7: usize = usize::MAX;
    let mut curr_i5: usize = usize::MAX;
    let mut i5_val: String;
    let mut delimiter = '\t';
    for line in lines {
        if line.len() < 5 {
            continue;
        }

        if header.len() == 0{
            header = line.to_lowercase().split('\t').map(|x| x.trim().to_string()).collect();
            if header.len() < 2{
                header = line.to_lowercase().split(',').map(|x| x.trim().to_string()).collect();
                delimiter = ',';
                if header.len() < 2{
                    panic!("Sample sheet columns should be separated by ',' or '\t'!");
                }
            }
            for header_itr in 0..header.len(){
                if header[header_itr] == "sample_id"{
                    curr_sample_id = header_itr;
                    sample_information.0[0] = header_itr;
                }else if header[header_itr] == "i7"{
                    curr_i7 = header_itr;
                    sample_information.0[1] = header_itr;
                }else if header[header_itr] == "i5"{
                    curr_i5 = header_itr;
                    sample_information.0[2] = header_itr;
                }else if header[header_itr] == "project_id"{
                    sample_information.0[3] = header_itr;
                }                      
            }
            if curr_sample_id == usize::MAX || curr_i7 == usize::MAX {
                panic!("sample_id and i7 must be in the sample sheet!");
            }
        }else{
            
            let vals: Vec<String> = line.split(delimiter).map(|x| x.trim().to_string()).collect();
            if vals[curr_i7] == "." || vals[curr_i7].len() < 3 {
                panic!("i7 ({}) should be longer than 3 chars!",  vals[curr_i7]);
            }
            
            if curr_i5 != usize::MAX {
                i5_val = vals[curr_i5].to_string();
                if i5_val == "." {
                    i5_val = String::new();
                }
                if i5_val != "" && i5_val.len() < 3 {
                    panic!(
                        "i5 ({}) should be longer than 3 chars!",
                        i5_val
                    );
                }

                sample_indexes.push([vals[curr_i7].to_string(), reverse_complement(&vals[curr_i7]).unwrap(), i5_val.to_string(), reverse_complement(&i5_val).unwrap()]);

            }else{
                sample_indexes.push([vals[curr_i7].to_string(), reverse_complement(&vals[curr_i7]).unwrap(), String::new(), String::new()]);
            }
            
            
            sample_information.1.push(vals.to_owned());
            
        }
        
    }

    Ok((sample_information.0, sample_information.1, sample_indexes))
}

pub fn add_match(curr_template: String, sample_itr: usize, matches_stat: &mut HashMap<String, Vec<usize>>, sample_list_ln: usize) {
    match matches_stat.get_mut(&curr_template){
        Some(sample_reads) => {
            sample_reads[sample_itr] += 1;
        },
        None => {
            let mut tmp = vec![0; sample_list_ln];
            tmp[sample_itr] += 1;
            matches_stat.insert(curr_template, tmp);
        }
    };
}

pub fn get_mgikit_template(initial_template:&String, umi: bool, barcode_length: usize, i7_len: usize, i5_len: usize) -> (String, String, String){
    
    let mut final_template: Vec<String> = Vec::new();
    let template_elements : Vec<String> = initial_template.split('_').map(String::from).collect();
    let rc_ls : Vec<String> = template_elements[1].split(':').map(String::from).collect();
    let mut index_ls : Vec<usize> = template_elements[0].split(':').map(|it| it.parse().unwrap()).collect();
    
    
        let mut shift = 0;
        let no_swap: bool = index_ls.len() == 1 || index_ls[0] < index_ls[1];
        index_ls.sort();

        let mut umi_loc = 0;
        let mut umi_size= 0;
        let mut item_added = 0;
        //println!("{:?}", index_ls);
        if index_ls[0] > 0{
            final_template.push(format!("--{}", index_ls[0]));
            shift += index_ls[0];
            umi_size = index_ls[0];
            item_added += 1;
        }
        //println!("1-   {:?}", final_template);
        if no_swap{
            final_template.push(format!("i7{}", i7_len));
            shift += i7_len;
        }else{
            final_template.push(format!("i5{}", i5_len));
            shift += i5_len;
        }
        item_added += 1;
        //println!("2-   {:?}", final_template);
        
        if index_ls.len() == 2{
            if shift < index_ls[1]{
                final_template.push(format!("--{}", index_ls[1] - shift));
                if index_ls[1] - shift > umi_size{
                    umi_loc = item_added;
                    umi_size = index_ls[1] - shift;
                }
                shift = index_ls[1];
                item_added += 1;
            }

            if no_swap{
                final_template.push(format!("i5{}", i5_len));
                shift += i5_len;
            }else{
                final_template.push(format!("i7{}", i7_len));
                shift += i7_len;
            }
            item_added += 1;
        }
        //println!("3-   {:?}", final_template);
        
        if shift < barcode_length{
            final_template.push(format!("--{}", barcode_length - shift));
            if barcode_length - shift > umi_size{
                umi_loc = item_added;
            }
        }
        //println!("4-   {:?}", final_template);
        
        if umi{
            //println!("{}  ->  {}", final_template[umi_loc], umi_loc);
            final_template[umi_loc] = final_template[umi_loc].replace("--", "um");
        }
        //println!("5-   {:?}", final_template);
        
    (final_template.join(":"), rc_ls[0].to_owned(), rc_ls[1].to_owned())
}

pub fn find_matches(sample_itr: usize, 
    matches_stat: &mut HashMap<String, Vec<usize>>, 
    sample_list_ln: usize, 
    read_barcode_seq: &String,
    sample_first_indx: &String,
    sample_first_indx_rc: &String,
    sample_second_indx: &String,
    sample_second_indx_rc: &String){
    
    //println!("{} - {} - {} - {}  -> {}", sample_first_indx, sample_first_indx_rc, sample_second_indx, sample_second_indx_rc, read_barcode_seq);
    let indexes  = match sample_second_indx.len() < 2 {
            true => vec![sample_first_indx, sample_first_indx_rc],
            false => vec![sample_first_indx, sample_first_indx_rc, sample_second_indx, sample_second_indx_rc]
    };
    //println!("{:?}", indexes);
        for i in 0..2{
            let first_matches: Vec<usize> = read_barcode_seq
            .match_indices(indexes[i])
            .map(|(i, _)| i)
            .collect();
            //println!("---------------------");
            //println!("{} -> {}, matches first: {:?}", i, indexes[i], first_matches);

            if first_matches.len() > 0{
                for j in 2..indexes.len(){
                    let second_matches: Vec<usize> = read_barcode_seq
                    .match_indices(indexes[j])
                    .map(|(i, _)| i )
                    .collect();
                   //println!("{} -> {}, matches second: {:?}", j, indexes[j], second_matches);
                    for second_match in &second_matches{
                        for first_match in &first_matches{
                            //println!("sample - {}: comparing  {}  ->  {}    {}:{}",sample_itr, first_match, second_match, i, j);

                            if *second_match >= first_match + sample_first_indx.len() || 
                                *first_match >= second_match + sample_second_indx.len(){
                                add_match(format!("{}:{}_{}:{}", 
                                                            first_match, 
                                                            second_match,
                                                            i,
                                                            j - 2), 
                                    sample_itr, 
                                    matches_stat, 
                                    sample_list_ln);
                            }
                            
                        }
                    }
                    
                    
                }
    
                if indexes.len() == 2{
                    for first_match in first_matches{
                        add_match(format!("{}_{}:.", 
                                                            first_match,
                                                            i), 
                                    sample_itr, 
                                    matches_stat, 
                                    sample_list_ln);
                    }   
                }
            }
        }
        
    
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use crate::samplesheet_utils::find_matches;

    use super::get_mgikit_template;
    
    #[test]
    fn test_get_mgikit_template() {
        let tests: [(&str, bool, usize, usize, usize); 10] = 
        [("0:8_0:0", true, 16, 8, 8), ("8_1:.", true, 16, 8, 0), ("0_1:.", false, 16, 8, 0), ("2_0:.", false, 10, 8, 0), 
        ("2:12_0:0", false, 20, 8, 8), ("2:10_0:0", false, 20, 8, 8), ("8:0_1:0", true, 24, 8, 8),
        ("0:16_0:0", true, 24, 8, 8), ("8:16_0:1", false, 24, 8, 8), ("8:16_0:1", true, 24, 8, 8)];
        let expected_out: [(&str, &str, &str); 10] = 
        [("i78:i58", "0", "0"), ("um8:i78", "1", "."), ("i78:--8", "1", "."), ("--2:i78", "0", "."), ("--2:i78:--2:i58", "0", "0"), 
        ("--2:i78:i58:--2", "0", "0"), ("i58:i78:um8", "1", "0"), ("i78:um8:i58", "0", "0"), ("--8:i78:i58", "0", "1"), 
        ("um8:i78:i58", "0", "1")];
        for i in 0..tests.len(){
            let res = get_mgikit_template(
                &String::from(tests[i].0), 
                tests[i].1,
                tests[i].2,
                tests[i].3,
                tests[i].4);

            assert_eq!(res.0, expected_out[i].0);
            assert_eq!(res.1, expected_out[i].1);
            assert_eq!(res.2, expected_out[i].2);
        }
        
    }
    
    #[test]
    fn test_find_matches() {
        let mut matches: HashMap<String, Vec<usize>> = HashMap::new();
        
        let tests:[(usize, usize, String, String, String, String, String, usize); 4] = [
            (0, 5, String::from("AAAAAAAACCCCCCCC"), String::from("CCCCCCCC"), String::from("GGGGGGGG"), String::from("TTTTTTTT"), String::from("AAAAAAAA"), 4), 
            (1, 5, String::from("CCACGGTCTGCGAGAG"), String::from("CCACGGTC"), String::from("GACCGTGG"), String::from("TGCGAGAG"), String::from("CTCTCGCA"), 5),
            (2, 5, String::from("ACGCTAGCTAGATTGA"), String::from("ACGCTAGCTA"), String::from("TAGCTAGCGT"), String::from(""), String::from(""), 6),
            (3, 7, String::from("CTCTACGCTAGCTAGATTGA"), String::from("ACGCTAGCTA"), String::from("TAGCTAGCGT"), String::from(""), String::from(""), 7)
            ];
        let expected: [(String, usize); 4] = [(String::from("8:0_0:1"), 0), (String::from("0:8_0:0"), 0), (String::from("0_0:."), 0), (String::from("4_0:."), 0)];
        
        for i in 0..tests.len(){
            for _ in 0..tests[i].7{
                find_matches(
                    tests[i].0, 
                    &mut matches,
                    tests[i].1,
                    &tests[i].2,
                    &tests[i].3,
                    &tests[i].4,
                    &tests[i].5,
                    &tests[i].6);
                    
            }
        }

        for i in 0..tests.len(){
            match matches.get(&expected[i].0){
                Some(template_info) => {
                    assert_eq!(template_info.len(), tests[i].1);
                    assert_eq!(template_info[tests[i].0], tests[i].7)
                },
                None => assert_eq!(0, 1)
            };
        }
        
    }
}