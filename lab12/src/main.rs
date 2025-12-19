use std::collections::HashMap;

fn main() {
    let motifs = vec![
        "GAGGTAAAC", 
        "TCCGTAAGT", 
        "CAGGTTGGA", 
        "ACAGTCAGT", 
        "TAGGTCATT", 
        "TAGGTACTG", 
        "ATGGTAACT", 
        "CAGGTATAC", 
        "TGTGTGAGT", 
        "AAGGTAAGT", 
    ];

    let sequence_s = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA";
    let motif_len = 9;
    let num_sequences = motifs.len() as f64;
    let bases = vec!['A', 'C', 'G', 'T'];

    let mut count_matrix: Vec<HashMap<char, u32>> = Vec::new();

    for i in 0..motif_len {
        let mut counts = HashMap::new();
        for &base in &bases {
            counts.insert(base, 0);
        }
        
        for motif in &motifs {
            let base = motif.chars().nth(i).unwrap();
            *counts.get_mut(&base).unwrap() += 1;
        }
        count_matrix.push(counts);
    }

    println!("1. Count Matrix:");
    print_matrix(&bases, &count_matrix, |val| format!("{:>4}", val));
    println!();

    let mut freq_matrix: Vec<HashMap<char, f64>> = Vec::new();

    for col in &count_matrix {
        let mut freqs = HashMap::new();
        for &base in &bases {
            let count = col[&base] as f64;
            let freq = count / num_sequences;
            freqs.insert(base, freq);
        }
        freq_matrix.push(freqs);
    }

    println!("3. Relative Frequencies Matrix (PPM):");
    print_matrix(&bases, &freq_matrix, |val| format!("{:>6.2}", val));
    println!();

    let null_model = 0.25;
    let mut log_matrix: Vec<HashMap<char, f64>> = Vec::new();

    for col in &freq_matrix {
        let mut logs = HashMap::new();
        for &base in &bases {
            let p_n = col[&base];
            if p_n > 0.0 {
                let score = (p_n / null_model).ln();
                logs.insert(base, score);
            } else {
                logs.insert(base, f64::NEG_INFINITY);
            }
        }
        log_matrix.push(logs);
    }

    println!("4. Log-Likelihoods Matrix (PWM):");
    print_matrix(&bases, &log_matrix, |val| {
        if val.is_infinite() {
            "  -inf".to_string()
        } else {
            format!("{:>6.2}", val)
        }
    });
    println!();

    println!("5. Analyze Sequence S: {}", sequence_s);
    println!("   Calculating scores for sliding window of size 9\n");

    let s_chars: Vec<char> = sequence_s.chars().collect();
  
    println!("{:<4} | {:<12} | {}", "Pos", "Window", "Score");
    println!("{:-<4}-+-{:-<12}-+-{:-<10}", "", "", "");

    let mut best_score = f64::NEG_INFINITY;
    let mut best_window = String::new();
    let mut best_pos = 0;

    for i in 0..=(s_chars.len() - motif_len) {
        let window = &s_chars[i..i+motif_len];
        let window_str: String = window.iter().collect();
        
        let mut current_score: f64 = 0.0;
        let mut impossible = false;

        for (pos, &base) in window.iter().enumerate() {
            let score = log_matrix[pos][&base];
            if score.is_infinite() {
                impossible = true;
                break;
            }
            current_score += score;
        }

        let score_display = if impossible {
            "-inf".to_string()
        } else {
            if current_score > best_score {
                best_score = current_score;
                best_window = window_str.clone();
                best_pos = i + 1;
            }
            format!("{:.4}", current_score)
        };

        println!("{:<4} | {} | {}", i + 1, window_str, score_display);
    }
    
    println!("\n---------------");
    println!("HIGHEST SIGNAL FOUND:");
    println!("Position: {}", best_pos);
    println!("Sequence: {}", best_window);
    println!("Score:    {:.4}", best_score);
    println!("----------------");
}

fn print_matrix<T, F>(bases: &Vec<char>, matrix: &Vec<HashMap<char, T>>, formatter: F)
where F: Fn(&T) -> String {
    print!("   ");
    for i in 1..=matrix.len() {
        print!("{:>6} ", i);
    }
    println!("\n-----------d-");
    for &base in bases {
        print!("{} |", base);
        for col in matrix {
            print!(" {} ", formatter(&col[&base]));
        }
        println!();
    }
}

// Do you have signals indicating that the S sequence contains an exon-intron border?
// The analysis identified a strong signal at Position 10 (AACGTAATC) with a score of 4.37. A positive score indicates that this sequence is significantly more likely to be an exon-intron boundary than a random DNA sequence (which would score 0 or negative).