const S: &str = "ABBA";

fn unique_kmers_in_order(s: &str, k: usize) -> Vec<String> {
    let n = s.len();
    if k == 0 || k > n {
        return Vec::new();
    }

    let bytes = s.as_bytes();
    let mut out: Vec<String> = Vec::new();

    let mut i: usize = 0;
    while i <= n - k {
        let window = &bytes[i..i + k];
        let kmer = std::str::from_utf8(window).unwrap();

        let mut seen = false;
        for existing in &out {
            if existing == kmer {
                seen = true;
                break;
            }
        }

        if !seen {
            out.push(kmer.to_string());
        }

        i += 1;
    }

    out
}

fn main() {
    let dinucs = unique_kmers_in_order(S, 2);
    let trinucs = unique_kmers_in_order(S, 3);

    println!("S=\"{}\"\n", S);

    println!("k=2 (dinucleotide):");
    if dinucs.is_empty() {
        println!(" ");
    } else {
        for km in &dinucs {
            println!("{}", km);
        }
    }
    println!();

    println!("k=3 (trinucleotide):");
    if trinucs.is_empty() {
        println!(" ");
    } else {
        for km in &trinucs {
            println!("{}", km);
        }
    }
}
