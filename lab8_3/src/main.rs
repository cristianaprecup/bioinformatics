use std::env;
use std::fs;

fn complement(base: char) -> char {
    match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        _ => 'N',
    }
}

fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement).collect()
}

fn read_fasta(path: &str) -> String {
    let content = fs::read_to_string(path).expect("Cannot read FASTA file");
    content
        .lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>()
}

fn find_inverted_repeats(seq: &str, min_len: usize, max_len: usize) {
    let chars: Vec<char> = seq.chars().collect();
    let n = chars.len();

    let max_spacer = 200;
    let mut count4 = 0;
    let mut count5 = 0;
    let mut count6 = 0;

    for len in min_len..=max_len {
        println!("Searching IR of length {}", len);

        for i in 0..=n - len {
            let left: String = chars[i..i + len].iter().collect();
            let rc = reverse_complement(&left);

            let end_j = (i + len + max_spacer).min(n - len);

            for j in i + len..=end_j {
                let right: String = chars[j..j + len].iter().collect();

                if right == rc {
                    println!(
                        "IR {} bp. {} at {} <-> {} at {}",
                        len, left, i, right, j
                    );

                    match len {
                        4 => count4 += 1,
                        5 => count5 += 1,
                        6 => count6 += 1,
                        _ => (),
                    }
                }
            }
        }
    }

    println!();
    println!("Summary:");
    println!("IR of length 4: {}", count4);
    println!("IR of length 5: {}", count5);
    println!("IR of length 6: {}", count6);
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        println!("Usage. cargo run <path_to_fasta>");
        return;
    }

    let fasta_path = &args[1];
    let seq = read_fasta(fasta_path);
    let seq_upper = seq.to_uppercase();

    println!("Loaded sequence with {} bases", seq_upper.len());
    println!("Searching for inverted repeats of length 4 to 6");

    find_inverted_repeats(&seq_upper, 4, 6);
}
