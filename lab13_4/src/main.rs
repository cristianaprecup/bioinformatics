use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Write};
use rand::Rng;
use serde::Deserialize;

#[derive(Deserialize)]
struct DnaAnalysis {
    #[allow(dead_code)] 
    sequence_length: usize,
    transition_matrix: HashMap<char, HashMap<char, f64>>,
}

#[derive(Deserialize)]
struct TextMarkovModel {
    #[allow(dead_code)]
    total_words_processed: usize,
    #[allow(dead_code)]
    vocabulary_size: usize,
    symbol_legend: HashMap<usize, String>, 
    transition_matrix: HashMap<usize, HashMap<usize, f64>>,
}

fn main() {
    println!("1. Generate DNA Sequence");
    println!("2. Generate English Text");
    print!("Select an engine (enter 1 or 2): ");
    io::stdout().flush().unwrap();

    let mut choice = String::new();
    io::stdin().read_line(&mut choice).expect("Failed to read input");

    match choice.trim() {
        "1" => run_dna_engine(),
        "2" => run_text_engine(),
        _ => println!("Invalid selection"),
    }
}

fn run_dna_engine() {    
    let file = match File::open("dna_transition_matrix.json") {
        Ok(f) => f,
        Err(_) => {
            eprintln!("Error: 'dna_transition_matrix.json' not found. Please run the DNA trainer first.");
            return;
        }
    };
    
    let model: DnaAnalysis = match serde_json::from_reader(file) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Error reading JSON: {}", e);
            return;
        }
    };

    let mut rng = rand::thread_rng();
    let bases: Vec<char> = model.transition_matrix.keys().cloned().collect();
    
    if bases.is_empty() { return; }

    let mut current_base = bases[rng.gen_range(0..bases.len())];
    let mut result = String::new();
    result.push(current_base);

    for _ in 0..49 {
        if let Some(transitions) = model.transition_matrix.get(&current_base) {
            let r: f64 = rng.r#gen();
            let mut cum_prob = 0.0;
            let mut next_base = current_base;

            for (target, &prob) in transitions {
                cum_prob += prob;
                if r <= cum_prob {
                    next_base = *target;
                    break;
                }
            }
            result.push(next_base);
            current_base = next_base;
        } else {
            break;
        }
    }

    println!("\nGenerated DNA:\n{}", result);
}

fn run_text_engine() {
    let file = match File::open("text_transitions.json") {
        Ok(f) => f,
        Err(_) => {
            eprintln!("Error: 'text_transitions.json' not found. Please run the Text trainer first.");
            return;
        }
    };

    let model: TextMarkovModel = match serde_json::from_reader(file) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("Error reading JSON: {}", e);
            return;
        }
    };

    let mut rng = rand::thread_rng();
    let all_ids: Vec<usize> = model.transition_matrix.keys().cloned().collect();
    
    if all_ids.is_empty() { return; }

    let mut current_id = all_ids[rng.gen_range(0..all_ids.len())];
    
    print!("\nGenerated Text:\n{}", model.symbol_legend.get(&current_id).unwrap());

    for _ in 0..100 {
        if let Some(transitions) = model.transition_matrix.get(&current_id) {
            let r: f64 = rng.r#gen();
            let mut cum_prob = 0.0;
            let mut next_id = current_id;
            let mut found = false;

            for (target, &prob) in transitions {
                cum_prob += prob;
                if r <= cum_prob {
                    next_id = *target;
                    found = true;
                    break;
                }
            }
            
            if !found && !transitions.is_empty() {
                next_id = *transitions.keys().last().unwrap();
            }

            if let Some(word) = model.symbol_legend.get(&next_id) {
                print!(" {}", word);
            }
            current_id = next_id;
        } else {
            current_id = all_ids[rng.gen_range(0..all_ids.len())];
            print!(". {}", model.symbol_legend.get(&current_id).unwrap());
        }
    }
    println!(".");
}