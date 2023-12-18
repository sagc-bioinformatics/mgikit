use std::collections::{HashMap, HashSet};
use itertools::Itertools;

/* 
pub fn get_mismatches(index: &String, allowed_mismatches: usize) -> HashSet<String> {
    let mut final_output = HashSet::new();
    let nec_ls = ['A', 'C', 'T', 'G', 'N'];
    let combination_ls = (0..index.len()).combinations(allowed_mismatches);

    for combination in combination_ls {
        let mut all_combinations: Vec<String> = nec_ls
            .iter()
            .filter(|&i| *i != index.chars().nth(combination[0]).unwrap())
            .map(|x| x.to_string())
            .collect();
        for nec_ind in 1..combination.len() {
            let ls_nec: Vec<char> = nec_ls
                .iter()
                .filter(|&i| *i != index.chars().nth(combination[nec_ind]).unwrap())
                .cloned()
                .collect();
            let mut tmp: Vec<String> = Vec::new();
            for x in all_combinations {
                for y in &ls_nec {
                    let mut z = x.clone();
                    z.push(*y);
                    tmp.push(z);
                }
            }
            all_combinations = tmp;
        }

        for combination_inner in all_combinations {
            let mut index_edited = index.clone();

            for nec_ind in 0..combination.len() {
                index_edited.replace_range(
                    combination[nec_ind]..(combination[nec_ind] + 1),
                    &combination_inner[nec_ind..(nec_ind + 1)],
                );
            }
            final_output.insert(index_edited.to_owned());
        }
    }

    final_output
}
*/
pub fn get_all_mismatches_old<'a>(
    all_indexes: &'a HashSet<String>,
    allowed_mismatches: usize,
) -> HashMap<String, (Vec<&'a String>, usize)> {
    if all_indexes.len() == 0{
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
        };

        for mismatch_itr in 1..(allowed_mismatches + 1) {
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
                                combination[nec_ind]..(combination[nec_ind] + 1),
                                &combination_inner[nec_ind..(nec_ind + 1)],
                            );
                        }
                        //println!("{}   {}", index, index_edited);

                        match final_output.get_mut(&index_edited) {
                            Some(x) => {
                                x[mismatch_itr].push(index);
                            }
                            None => {
                                let mut tmp = Vec::new();
                                for mis_itr in 0..(allowed_mismatches + 1) {
                                    if mis_itr == mismatch_itr {
                                        tmp.push(vec![index]);
                                    } else {
                                        tmp.push(Vec::new());
                                    }
                                }
                                final_output.insert(index_edited.clone(), tmp);
                            }
                        };
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
        for i in 0..(allowed_mismatches + 1) {
            if value[i].len() > 0 {
                //let mut tmp_z = value[i].clone();
                //tmp_z.sort_by(|a, b| a.cmp(b));
                sample_index_map.insert(key, (
                    value[i].to_owned(), 
                    i)
                );
                break;
            }
        }
    }

    sample_index_map
}


pub fn get_all_mismatches<'a>(
    all_indexes: &'a HashSet<String>,
    allowed_mismatches: usize,
) -> HashMap<Vec<u8>, (Vec<&'a String>, usize)> {
    if all_indexes.len() == 0{
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
        };

        for mismatch_itr in 1..(allowed_mismatches + 1) {
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
                                combination[nec_ind]..(combination[nec_ind] + 1),
                                &combination_inner[nec_ind..(nec_ind + 1)],
                            );
                        }
                        //println!("{}   {}", index, index_edited);

                        match final_output.get_mut(&index_edited) {
                            Some(x) => {
                                x[mismatch_itr].push(index);
                            }
                            None => {
                                let mut tmp = Vec::new();
                                for mis_itr in 0..(allowed_mismatches + 1) {
                                    if mis_itr == mismatch_itr {
                                        tmp.push(vec![index]);
                                    } else {
                                        tmp.push(Vec::new());
                                    }
                                }
                                final_output.insert(index_edited.clone(), tmp);
                            }
                        };
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
        for i in 0..(allowed_mismatches + 1) {
            if value[i].len() > 0 {
                //let mut tmp_z = value[i].clone();
                //tmp_z.sort_by(|a, b| a.cmp(b));
                sample_index_map.insert(key.as_bytes().to_vec(), (
                    value[i].to_owned(), 
                    i)
                );
                break;
            }
        }
    }

    sample_index_map
}
