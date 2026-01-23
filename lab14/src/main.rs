use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Base {
    A, C, G, T,
}

impl Base {
    fn from_char(c: char) -> Option<Self> {
        match c {
            'A' => Some(Base::A),
            'C' => Some(Base::C),
            'G' => Some(Base::G),
            'T' => Some(Base::T),
            _ => None,
        }
    }

    fn index(&self) -> usize {
        match self {
            Base::A => 0,
            Base::C => 1,
            Base::G => 2,
            Base::T => 3,
        }
    }
    
    fn all() -> [Base; 4] {
        [Base::A, Base::C, Base::G, Base::T]
    }
}

type Matrix4x4 = [[f64; 4]; 4];

fn main() {
    let s1_pos = "ATCGATTCGATATCATACACGTAT"; 
    let s2_neg = "CTCGACTAGTATGAAGTCCACGCTTG"; 
    let s_target = "CAGGTTGGAAACGTAA"; 

    let model_plus = train_markov_model(s1_pos);
    let model_minus = train_markov_model(s2_neg);
    
    println!("(+) Model Matrix Frequency:");
    print_matrix(&model_plus);
    
    println!("\n(-) Model Matrix Frequency:");
    print_matrix(&model_minus);

    let llr_matrix = calculate_log_likelihood(&model_plus, &model_minus);
    print_matrix(&llr_matrix);

    let score = score_sequence(s_target, &llr_matrix);
    
    println!("Sequence S: {}", s_target);
    println!("Total score: {:.4}", score);
    
    if score > 0.0 {
        println!("Conclusion: Sequence S BELONGS to a CpG island.");
    } else {
        println!("Conclusion: Sequence S does NOT belong to a CpG island.");
    }
}

fn train_markov_model(sequence: &str) -> Matrix4x4 {
    let bases: Vec<Base> = sequence.chars().filter_map(Base::from_char).collect();
    let mut counts = [[0.0; 4]; 4];
    let mut row_sums = [0.0; 4];

    for window in bases.windows(2) {
        let from = window[0].index();
        let to = window[1].index();
        counts[from][to] += 1.0;
        row_sums[from] += 1.0;
    }

    let mut probabilities = [[0.0; 4]; 4];
    for r in 0..4 {
        for c in 0..4 {
            if row_sums[r] > 0.0 {
                probabilities[r][c] = counts[r][c] / row_sums[r];
            } else {
                probabilities[r][c] = 0.0;
            }
        }
    }
    probabilities
}

fn calculate_log_likelihood(plus: &Matrix4x4, minus: &Matrix4x4) -> Matrix4x4 {
    let mut llr = [[0.0; 4]; 4];

    for r in 0..4 {
        for c in 0..4 {
            let p_plus = plus[r][c];
            let p_minus = minus[r][c];

            if p_plus == 0.0 {
                llr[r][c] = 0.0; 
            } else if p_minus == 0.0 {
                llr[r][c] = 0.0; 
            } else {
                llr[r][c] = (p_plus / p_minus).log2();
            }
        }
    }
    llr
}

fn score_sequence(seq: &str, llr_matrix: &Matrix4x4) -> f64 {
    let bases: Vec<Base> = seq.chars().filter_map(Base::from_char).collect();
    let mut total_score = 0.0;

    print!("Transitions: ");
    for window in bases.windows(2) {
        let from = window[0];
        let to = window[1];
        let score = llr_matrix[from.index()][to.index()];
        total_score += score;
    }
    total_score
}

fn print_matrix(matrix: &Matrix4x4) {
    println!("\tA\tC\tG\tT");
    let rows = Base::all();
    for (i, row_label) in rows.iter().enumerate() {
        print!("{:?}\t", row_label);
        for val in matrix[i].iter() {
            print!("{:.3}\t", val);
        }
        println!();
    }
}