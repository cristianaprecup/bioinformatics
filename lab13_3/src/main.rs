use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use serde::Serialize;
use regex::Regex;

#[derive(Serialize)]
struct TextMarkovModel {
    total_words_processed: usize,
    vocabulary_size: usize,
    symbol_legend: HashMap<usize, String>, 
    transition_matrix: HashMap<usize, HashMap<usize, f64>>,
}

fn main() -> io::Result<()> {
    let path = "english_text.txt";
    let file = File::open(path).expect("Could not open 'english_text.txt'.");
    let reader = BufReader::new(file);

    let re = Regex::new(r"[^a-zA-Z0-9\s]").unwrap();
    let mut raw_words: Vec<String> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() { continue; }
        
        let cleaned = re.replace_all(&line, "").to_lowercase();
        for word in cleaned.split_whitespace() {
            raw_words.push(word.to_string());
        }
    }

    if raw_words.len() < 2 {
        eprintln!("Error: Text is too short to calculate transitions.");
        return Ok(());
    }

    println!("Processed {} words.", raw_words.len());

    let mut word_to_id: HashMap<String, usize> = HashMap::new();
    let mut id_to_word: HashMap<usize, String> = HashMap::new();
    let mut next_id: usize = 0;

    let mut sequence_ids: Vec<usize> = Vec::with_capacity(raw_words.len());

    for word in raw_words {
        let id = if let Some(&id) = word_to_id.get(&word) {
            id
        } else {
            let id = next_id;
            word_to_id.insert(word.clone(), id);
            id_to_word.insert(id, word);
            next_id += 1;
            id
        };
        sequence_ids.push(id);
    }

    let mut counts: HashMap<usize, HashMap<usize, f64>> = HashMap::new();

    for window in sequence_ids.windows(2) {
        let current_id = window[0];
        let next_id = window[1];

        counts.entry(current_id)
            .or_insert_with(HashMap::new)
            .entry(next_id)
            .and_modify(|c| *c += 1.0)
            .or_insert(1.0);
    }

    let mut probability_matrix: HashMap<usize, HashMap<usize, f64>> = HashMap::new();

    for (source_id, targets) in counts {
        let total_transitions: f64 = targets.values().sum();
        let mut row = HashMap::new();

        for (target_id, count) in targets {
            let prob = count / total_transitions;
            row.insert(target_id, (prob * 10000.0).round() / 10000.0);
        }
        probability_matrix.insert(source_id, row);
    }

    let model = TextMarkovModel {
        total_words_processed: sequence_ids.len(),
        vocabulary_size: id_to_word.len(),
        symbol_legend: id_to_word,
        transition_matrix: probability_matrix,
    };

    let output_file = "text_transitions.json";
    let file = File::create(output_file)?;
    serde_json::to_writer_pretty(file, &model)?;

    println!("Analysis saved to '{}'", output_file);
    println!("Vocabulary Size: {} unique words", model.vocabulary_size);

    Ok(())
}