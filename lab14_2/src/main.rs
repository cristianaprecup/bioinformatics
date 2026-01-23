use std::collections::{HashMap, HashSet};
use std::f64;

const EMINESCU_TEXT: &str = "
Somnoroasele păsări
Pe la cuiburi se adună,
Se ascund în ramurele -
Noapte bună!

Doar izvoarele suspină,
Pe când codrul negru tace;
Dorm și florile-n grădină -
Dormi în pace!

Trece lebăda pe ape
Între trestii să se culce -
Fie-ți îngerii aproape,
Somnul dulce!

Peste-a nopții feerie
Se ridică mândra lună,
Totu-i vis și armonie -
Noapte bună!
";

const STANESCU_TEXT: &str = "
Leoaică tânără, iubirea
mi-ai sărit în față.
Mă pândise-n încordare
mai demult.
Colții albi mi i-a înfipt în față,
m-a mușcat leoaica, azi, de față.
Și deodată-n jurul meu, natura
se făcu un cerc, de-a-dura,
când mai larg, când mai aproape,
ca o strângere de ape.
Și privirea-n sus țâșni,
curcubeu tăiat în două,
și auzul o-ntâlni
tocmai lângă ciocârlii.

Mi-am dus mâna la sprânceană,
la tâmplă și la bărbie,
dar mâna nu le mai știe.
Și alunecă în neștire
pe-un deșert în strălucire,
peste care trece-alene
o leoaică arămie
cu mișcările viclene,
încă-o vreme,
și-ncă-o vreme...
";

const TILT_TEXT: &str = "
Somnoroasele păsări se adună pe la cuiburi.
Noapte bună, dormi în pace!
Leoaică tânără cu mișcările viclene
mi-ai sărit în față încă-o vreme.
";

fn main() {
    let mut vocab: HashSet<String> = HashSet::new();
    let eminescu_words = tokenize(EMINESCU_TEXT);
    let stanescu_words = tokenize(STANESCU_TEXT);
    
    for w in &eminescu_words { vocab.insert(w.clone()); }
    for w in &stanescu_words { vocab.insert(w.clone()); }

    let mut word_to_char: HashMap<String, char> = HashMap::new();
    let mut char_to_word: HashMap<char, String> = HashMap::new();
    
    for (i, word) in vocab.iter().enumerate() {
        let symbol = std::char::from_u32((i + 33) as u32).unwrap_or('?');
        word_to_char.insert(word.clone(), symbol);
        char_to_word.insert(symbol, word.clone());
    }

    println!("Vocabulary Size: {} unique words encoded as characters.", vocab.len());

    let seq_eminescu = encode_sequence(&eminescu_words, &word_to_char);
    let seq_stanescu = encode_sequence(&stanescu_words, &word_to_char);

    let all_chars: Vec<char> = word_to_char.values().cloned().collect();
    
    let prob_matrix_eminescu = train_markov(&seq_eminescu, &all_chars);
    let prob_matrix_stanescu = train_markov(&seq_stanescu, &all_chars);

    let llm = build_log_likelihood_matrix(&prob_matrix_eminescu, &prob_matrix_stanescu, &all_chars);

    println!("Analyzing Tilt Text:");

    // ... inside main() ...

    // 5. Scan the Tilt Text
    let tilt_words = tokenize(TILT_TEXT);
    let seq_tilt_str = encode_sequence(&tilt_words, &word_to_char);
    
    // CRITICAL FIX: Convert String to Vec<char> to handle multi-byte chars safely
    let seq_tilt: Vec<char> = seq_tilt_str.chars().collect();
    
    let window_size = 4;
    let mut scores = Vec::new();
    let mut x_labels = Vec::new();

    if seq_tilt.len() < window_size {
        println!("Text too short for analysis.");
        return;
    }

    for window in seq_tilt.windows(window_size) {
        let score = score_window(window, &llm);
        scores.push(score);
        
        let first_char = window[0];
        let start_word = char_to_word.get(&first_char).unwrap();
        x_labels.push(start_word.clone());
    }

    draw_ascii_chart(&scores, &x_labels);
}


fn tokenize(text: &str) -> Vec<String> {
    text.to_lowercase()
        .chars()
        .map(|c| if c.is_alphanumeric() || c == '-' { c } else { ' ' })
        .collect::<String>()
        .split_whitespace()
        .map(|s| s.to_string())
        .collect()
}

fn encode_sequence(words: &[String], map: &HashMap<String, char>) -> String {
    words.iter()
        .map(|w| *map.get(w).unwrap_or(&'?'))
        .collect()
}

fn train_markov(seq: &str, alphabet: &[char]) -> HashMap<(char, char), f64> {
    let mut counts: HashMap<(char, char), f64> = HashMap::new();
    let mut row_sums: HashMap<char, f64> = HashMap::new();
    let chars: Vec<char> = seq.chars().collect();

    let alpha = 0.1; 
    let vocab_size = alphabet.len() as f64;

    for &c1 in alphabet {
        row_sums.insert(c1, alpha * vocab_size); 
        for &c2 in alphabet {
            counts.insert((c1, c2), alpha);
        }
    }

    for window in chars.windows(2) {
        let from = window[0];
        let to = window[1];
        
        if from == '?' || to == '?' { continue; }

        *counts.get_mut(&(from, to)).unwrap() += 1.0;
        *row_sums.get_mut(&from).unwrap() += 1.0;
    }

    let mut probs = HashMap::new();
    for ((from, to), count) in counts {
        let total = row_sums.get(&from).unwrap();
        probs.insert((from, to), count / total);
    }
    probs
}

fn build_log_likelihood_matrix(
    model_pos: &HashMap<(char, char), f64>, 
    model_neg: &HashMap<(char, char), f64>,
    alphabet: &[char]
) -> HashMap<(char, char), f64> {
    let mut llm = HashMap::new();
    
    for &from in alphabet {
        for &to in alphabet {
            let p_pos = model_pos.get(&(from, to)).unwrap_or(&1e-10); // Eminescu
            let p_neg = model_neg.get(&(from, to)).unwrap_or(&1e-10); // Stanescu
            
            let score = p_pos.ln() - p_neg.ln();
            llm.insert((from, to), score);
        }
    }
    llm
}

fn score_window(window: &[char], llm: &HashMap<(char, char), f64>) -> f64 {
    let mut score = 0.0;
    
    for w in window.windows(2) {
        let from = w[0];
        let to = w[1];
        score += llm.get(&(from, to)).unwrap_or(&0.0);
    }
    score
}


fn draw_ascii_chart(data: &[f64], labels: &[String]) {
    println!("\n--- Authorship Trend Chart ---");
    println!("(+) Positive = Eminescu | (-) Negative = Stanescu\n");

    let max_val = data.iter().cloned().fold(f64::MIN, f64::max).max(1.0);
    let min_val = data.iter().cloned().fold(f64::MAX, f64::min).min(-1.0);
    let range = max_val - min_val;
    let height = 10;

    for r in 0..=height {
        let current_y = max_val - (r as f64 * range / height as f64);
        print!("{:>5.1} | ", current_y);

        for &val in data {
            let threshold = range / (2.0 * height as f64);
            if (val - current_y).abs() < threshold {
                if val > 0.0 { print!("+"); } else { print!("-"); }
            } else if (val > current_y && r == height) || (val < current_y && r == 0) {
                 print!(" ");
            } else {
                print!(" ");
            }
            print!("   "); 
        }
        println!();
    }
    
    print!("      ");
    for _ in 0..data.len() { print!("----"); }
    println!();
    
    print!("      ");
    for label in labels {
        let short = if label.len() > 3 { &label[0..3] } else { label };
        print!("{:<3} ", short);
    }
    println!("\n");
    
    let avg: f64 = data.iter().sum::<f64>() / data.len() as f64;
    println!("Average Score: {:.4}", avg);
    if avg > 0.5 { println!("Result: Text leans towards MIHAI EMINESCU."); }
    else if avg < -0.5 { println!("Result: Text leans towards NICHITA STANESCU."); }
    else { println!("Result: Hybrid or Ambiguous text."); }
}