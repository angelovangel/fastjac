
use clap::{App, Arg};
use std::{collections::{HashSet}, hash::Hash};
extern crate clap;
//use std::cmp::Ordering;
use kseq::parse_path;


fn jaccard<K>(h1: HashSet<K>, h2: HashSet<K>) -> Vec<f32> where K: Hash+Eq{
    //let set1 = h1.keys().
    let i = h1.intersection(&h2).count() as f32;
    let u = h1.union(&h2).count() as f32;
    let l1 = h1.len() as f32;
    let l2 = h2.len() as f32;


    let jac = i / u; // jaccard
    //let olap = i / s as f32; // overlap
    let sd = (2.0 * i) / (h1.len() as f32 + h2.len() as f32); //Sorensen Dice
    let cont = i / h1.len() as f32; // containment
    return vec![l1, l2, i, u, jac, sd, cont];
}

fn main() {
    let argmatches = App::new("fastjac")
    .version("0.1.0")
    .author("https://github.com/angelovangel")
    .about("get various kmer distance metrics (Jaccard, Overlap, Sørensen-Dice) between fastq/fasta files")

    .arg(Arg::with_name("kmer")
    .long("kmer")
    .short("k")
    .required(true)
    .takes_value(true)
    .help("k-mer size to use"))
    
    .arg(Arg::with_name("query")
    .long("query")
    .short("q")
    .takes_value(true)
    .required(true)
    .help("Path to query fastx file"))
    
    .arg(Arg::with_name("valid")
    .long("valid")
    .short("v")
    .required(false)
    .takes_value(false)
    .help("Remove k-mers with ambigous bases, e.g. only the ones containing ATGC are retained"))
    //.index(1))
    
    .arg(Arg::with_name("ref")
    .long("ref")
    .short("r")
    .takes_value(true)
    .required(true)
    .help("Path to reference fastx file"))

    //.index(2))
    .get_matches();
    
    //

    let k = argmatches.value_of("kmer").unwrap().trim().parse::<usize>().expect("k-mer length argument not valid!");
    //let mut kmer_counts: HashMap<String,i32> = HashMap::new();

    let mut kmer_seq: HashSet<String> = HashSet::new();
    let mut kmer_ref: HashSet<String> = HashSet::new();
 
    // process ref
    let file2 = argmatches.value_of("ref").unwrap().to_string();
    let mut records_r = parse_path(file2).unwrap();

    while let Some(record) = records_r.iter_record().unwrap() {

        let ref_str = record.seq();
        for c in 0..ref_str.len() - k + 1 {
            let subseq = &ref_str[c..c + k];
            if argmatches.is_present("valid") {
                if subseq.chars().all(|x| matches!(x, 'A'|'T'|'G'|'C'|'a'|'t'|'g'|'c') ) {
                    kmer_ref.insert( subseq.to_string() );
                }
            } else {
                kmer_ref.insert(subseq.to_string());
            }
        }
    }

    // process query
    let file1 = argmatches.value_of("query").unwrap().to_string();
    let mut records_q = parse_path(file1).unwrap();    

    while let Some(record) = records_q.iter_record().unwrap() {

        let query_str = record.seq();
        if query_str.len() < k {
                continue; // skip if query is shorter than k
            }
        for c in 0..query_str.len() - k + 1 {
            let subseq = &query_str[c..c + k];
        
            if argmatches.is_present("valid") {  
                if subseq.chars().all(|x| matches!(x, 'A'|'T'|'G'|'C'|'a'|'t'|'g'|'c') ) {
                    kmer_seq.insert( subseq.to_string() );
                }
            } else {
                kmer_seq.insert( subseq.to_string() );
            }
        }
    }

    //println!("{:?}", kmer_seq );
    let distance = jaccard(kmer_seq, kmer_ref);
    
    println!("len_query\tlen_ref\tintersection\tunion\tjaccard\tSD\tcontainment");
    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", distance[0], distance[1], distance[2], distance[3], distance[4], distance[5], distance[6]);

}
