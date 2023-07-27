
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


/*
fn hamming_distance(sub_seq:&[u8], seq:&[u8], start:usize, length:usize) -> u32{
    let mut mismatches: u32 = 0;
    for i in 0..length{
        //println!("{} vs {}  ->  {}   {}", i, i+start, sub_seq[i], seq[i + start]);
        if sub_seq[i] != seq[i + start]{
            mismatches += 1;
        }
    }
    mismatches
}
*/


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