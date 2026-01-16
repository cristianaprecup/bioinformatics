use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use serde::Serialize;

const DNA_BASES: [char; 4] = ['A', 'C', 'G', 'T'];
/// Struct to represent the final JSON output
#[derive(Serialize)]
struct DnaAnalysis {
    sequence_length: usize,
    transition_matrix: HashMap<char, HashMap<char, f64>>,
}

fn main() {
    let dna_sequence = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAA"; 
    
    println!("Analyzing Sequence: {}", dna_sequence);

    let mut counts: HashMap<char, HashMap<char, f64>> = HashMap::new();
    
    for &source in &DNA_BASES {
        let mut row = HashMap::new();
        for &target in &DNA_BASES {
            row.insert(target, 0.0);
        }
        counts.insert(source, row);
    }

    let chars: Vec<char> = dna_sequence.chars().collect();
    
    if chars.len() < 2 {
        eprintln!("Sequence too short to calculate transitions.");
        return;
    }

    for window in chars.windows(2) {
        let current_base = window[0];
        let next_base = window[1];

        if isValidBase(current_base) && isValidBase(next_base) {
             if let Some(row) = counts.get_mut(&current_base) {
                if let Some(count) = row.get_mut(&next_base) {
                    *count += 1.0;
                }
            }
        }
    }

    let mut probability_matrix: HashMap<char, HashMap<char, f64>> = HashMap::new();

    for &source in &DNA_BASES {
        let row_counts = counts.get(&source).unwrap();
        let total_occurrences: f64 = row_counts.values().sum();
        
        let mut probability_row = HashMap::new();

        for &target in &DNA_BASES {
            let count = row_counts.get(&target).unwrap();
            let probability = if total_occurrences > 0.0 {
                count / total_occurrences
            } else {
                0.0
            };
            
            probability_row.insert(target, (probability * 1000.0).round() / 1000.0);
        }
        probability_matrix.insert(source, probability_row);
    }

    let output_data = DnaAnalysis {
        sequence_length: chars.len(),
        transition_matrix: probability_matrix,
    };

    let filename = "dna_transition_matrix.json";
    match save_to_json(&output_data, filename) {
        Ok(_) => println!("Matrix saved to '{}'", filename),
        Err(e) => eprintln!("Error saving file: {}", e),
    }
}

fn isValidBase(c: char) -> bool {
    matches!(c, 'A' | 'C' | 'G' | 'T')
}

fn save_to_json(data: &DnaAnalysis, filename: &str) -> std::io::Result<()> {
    let file = File::create(filename)?;
    serde_json::to_writer_pretty(file, data)?;
    Ok(())
}