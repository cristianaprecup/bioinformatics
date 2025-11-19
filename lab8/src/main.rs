use rand::{Rng, SeedableRng};
use rand::rngs::StdRng;

const BASE_DNA_LEN: usize = 250;

const TRANSPOSONS: [&str; 4] = [
    "ATGCGTACGA",   // T1
    "TTACGTTACG",   // T2
    "CGTACGCGTA",   // T3
    "GATTACAGAT",   // T4
];

fn main() {
    let base_dna = generate_random_dna(BASE_DNA_LEN);
    println!("Initial sequence, length {}:", base_dna.len());
    println!("{base_dna}\n");

    let (dna_with_tps, inserted_positions) = insert_transposons(base_dna, &TRANSPOSONS);

    println!(
        "Final sequence, length {}:",
        dna_with_tps.len()
    );
    println!("{dna_with_tps}\n");

    println!("Real positions of inserted transposons:");
    for (name, start, end) in &inserted_positions {
        println!("{name}: start = {start}, end = {end}");
    }
    println!();

    let detected = detect_transposons(&dna_with_tps, &TRANSPOSONS);

    println!("Detected positions by the algorithm:");
    for (pattern, start, end) in detected {
        println!("{pattern}: start = {start}, end = {end}");
    }
}

fn generate_random_dna(len: usize) -> String {
    let mut rng = StdRng::seed_from_u64(42); 
    let bases = ['A', 'C', 'G', 'T'];

    (0..len)
        .map(|_| {
            let idx = rng.gen_range(0..bases.len());
            bases[idx]
        })
        .collect()
}

fn insert_transposons(
    base_dna: String,
    transposons: &[&str],
) -> (String, Vec<(String, usize, usize)>) {
    let mut seq: Vec<char> = base_dna.chars().collect();
    let mut positions: Vec<(String, usize, usize)> = Vec::new();

    let planned_positions: [usize; 4] = [
        50,   // T1
        120,  // T2
        55,   // T3
        180,  // T4
    ];

    let mut offset = 0;

    for (i, &pattern) in transposons.iter().enumerate() {
        let original_pos = planned_positions[i];

        let real_pos = original_pos + offset;

        let pattern_chars: Vec<char> = pattern.chars().collect();
        let len_p = pattern_chars.len();

        seq.splice(real_pos..real_pos, pattern_chars.clone());

        let start = real_pos;
        let end = real_pos + len_p - 1;

        positions.push((format!("T{}", i + 1), start, end));

        offset += len_p;
    }

    let final_seq: String = seq.iter().collect();
    (final_seq, positions)
}

fn detect_transposons(
    dna: &str,
    transposons: &[&str],
) -> Vec<(String, usize, usize)> {
    let mut results = Vec::new();

    for (i, &pattern) in transposons.iter().enumerate() {
        for (start, _) in dna.match_indices(pattern) {
            let end = start + pattern.len() - 1;
            results.push((format!("T{} ({})", i + 1, pattern), start, end));
        }
    }

    results
}
