fn read_dna_from_stdin() -> Vec<u8> {
    let mut input = String::new();
    io::stdin()
        .read_to_string(&mut input)
        .expect("cannot read from stdin");

    let mut seq = String::new();
    for line in input.lines() {
        if line.starts_with('>') {
            continue;
        }
        for ch in line.chars() {
            let up = ch.to_ascii_uppercase();
            if matches!(up, 'A' | 'C' | 'G' | 'T') {
                seq.push(up);
            }
        }
    }

    seq.into_bytes()
}

fn find_cuts(dna: &[u8], enzyme: &Enzyme) -> Vec<usize> {
    let site_bytes = enzyme.site.as_bytes();
    let site_len = site_bytes.len();
    let mut cuts = Vec::new();

    if site_len == 0 || site_len > dna.len() {
        return cuts;
    }

    let mut i = 0;
    while i + site_len <= dna.len() {
        if &dna[i..i + site_len] == site_bytes {
            let cut_pos = i + enzyme.cut_offset;
            cuts.push(cut_pos);
        }
        i += 1;
    }

    cuts.sort_unstable();
    cuts.dedup();
    cuts
}

fn compute_fragments(seq_len: usize, cuts: &[usize]) -> Vec<usize> {
    let mut boundaries = Vec::new();
    boundaries.push(0);
    boundaries.extend_from_slice(cuts);
    boundaries.push(seq_len);

    boundaries.sort_unstable();
    boundaries.dedup();

    let mut fragments = Vec::new();
    for w in boundaries.windows(2) {
        let start = w[0];
        let end = w[1];
        if end > start {
            fragments.push(end - start);
        }
    }

    fragments
}

fn simulate_gel(results: &[(String, Vec<usize>)]) {
    if results.is_empty() {
        println!("Nu există rezultate pentru gel.");
        return;
    }

    let largest_lane_height = 30;

    let mut min_len = usize::MAX;
    let mut max_len = 0;

    for (_, frags) in results {
        for &f in frags {
            if f < min_len { min_len = f; }
            if f > max_len { max_len = f; }
        }
    }

    let lane_width = 12;

    let mut matrix = vec![vec![' '; lane_width * results.len()]; largest_lane_height];

    for (lane_index, (_, frags)) in results.iter().enumerate() {

        for &f in frags {
            let relative = (f - min_len) as f64 / (max_len - min_len) as f64;
            let y = (relative * (largest_lane_height as f64 - 2.0)) as usize;

            let row = largest_lane_height - 2 - y;
            let col_start = lane_index * lane_width + 2;

            for x in 0..6 {
                if col_start + x < matrix[row].len() {
                    matrix[row][col_start + x] = '#';
                }
            }
        }
    }

    println!("\nGel of DNA fragments (ASCII simulation)\n");
    println!("      ↓ DNA migration");
    println!("      +");

    for row in &matrix {
        print!("|");
        for ch in row {
            print!("{}", ch);
        }
        println!("|");
    }

    print!("|");
    for _ in 0..(lane_width * results.len()) {
        print!("=");
    }
    println!("|");

    for (i, (name, _)) in results.iter().enumerate() {
        let pad = lane_width * i + 2;
        for _ in 0..pad { print!(" "); }
        println!("{}", name);
    }

    println!();
}


use std::env;
use std::fs;
use std::io::{self, Read};

#[derive(Debug)]
struct Enzyme {
    name: &'static str,
    site: &'static str,
    cut_offset: usize,
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let fasta_path = &args[1];

    let fasta_content = fs::read_to_string(fasta_path)
        .expect("cannot read FASTA file");

    let dna = parse_fasta(&fasta_content);
    let len = dna.len();

    println!("Sequence length: {} nucleotides\n", len);

    let enzymes = vec![
        Enzyme { name: "EcoRI",  site: "GAATTC", cut_offset: 1 },
        Enzyme { name: "BamHI",  site: "GGATCC", cut_offset: 1 },
        Enzyme { name: "HindIII",site: "AAGCTT", cut_offset: 1 },
        Enzyme { name: "TaqI",   site: "TCGA",   cut_offset: 1 },
        Enzyme { name: "HaeIII", site: "GGCC",   cut_offset: 2 },
    ];

    let mut all_results = Vec::new();

    for enzyme in &enzymes {
        let cuts = find_cuts(&dna, enzyme);
        let fragments = compute_fragments(len, &cuts);
        println!("=== {} ===", enzyme.name);
        println!("Site: {}", enzyme.site);
        println!("Cuts: {}", cuts.len());
        println!("Positions: {:?}", cuts.iter().map(|x| x + 1).collect::<Vec<_>>());
        println!("Fragments: {:?}\n", fragments);
        all_results.push((enzyme.name.to_string(), fragments));
    }

    println!("=== Simulate gel ===\n");
    simulate_gel(&all_results);
}

fn parse_fasta(content: &str) -> Vec<u8> {
    let mut seq = Vec::new();
    for line in content.lines() {
        if line.starts_with('>') { continue; }
        for ch in line.chars() {
            let up = ch.to_ascii_uppercase();
            if matches!(up, 'A' | 'C' | 'G' | 'T') {
                seq.push(up as u8);
            }
        }
    }
    seq
}
