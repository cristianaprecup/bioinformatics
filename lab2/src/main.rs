const SEQ: &str = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA";

fn generate_all_k(k: usize) -> Vec<String> {
    let mut result = vec![String::new()];
    for _ in 0..k {
        let mut next = Vec::with_capacity(result.len() * 4);
        for prefix in result {
            for &c in &['A', 'C', 'G', 'T'] {
                let mut s = prefix.clone();
                s.push(c);
                next.push(s);
            }
        }
        result = next;
    }
    result
}

fn count_overlapping(text: &str, pattern: &str) -> usize {
    if pattern.is_empty() {
        return 0;
    }
    let tb = text.as_bytes();
    let pb = pattern.as_bytes();
    let n = tb.len();
    let m = pb.len();
    if m > n {
        return 0;
    }
    let mut count = 0;
    for i in 0..=(n - m) {
        let window = &tb[i..i + m];
        if window == pb {
            count += 1;
        }
    }
    count
}

fn compute_per_k(k: usize) -> (Vec<(String, usize, usize, f64)>, usize) {
    let n = SEQ.len();
    let total_windows: usize;
    if k <= n {
        total_windows = n - k + 1;
    } else {
        total_windows = 0;
    }
    let mut rows: Vec<(String, usize, usize, f64)> = Vec::new();
    let mut kmers = generate_all_k(k);
    kmers.sort();
    for kmer in kmers {
        let c = count_overlapping(SEQ, &kmer);
        let pct = if total_windows > 0 {
            (c as f64) * 100.0 / (total_windows as f64)
        } else {
            0.0
        };
        rows.push((kmer, c, total_windows, pct));
    }
    (rows, total_windows)
}

fn print_k_percentages(k: usize) {
    let (rows, total) = compute_per_k(k);
    println!("combination, count, total, percentage");
    for (kmer, count, _total, pct) in rows {
        println!("{}, {}, {}, {:.3}", kmer, count, total, pct);
    }
    println!();
}

fn main() {
    println!("for dinucleotides:");
    print_k_percentages(2);
    println!("for trinucleotides:");
    print_k_percentages(3);
}
