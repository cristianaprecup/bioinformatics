use std::{env, fs, process::exit};

const W: usize = 8;         
const NA_MOLAR: f64 = 0.05;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!("usage: {} <input.fasta>", args[0]);
        exit(1);
    }

    let content = fs::read_to_string(&args[1]).expect("cannot read file");
    let entries = parse_fasta(&content);
    if entries.is_empty() {
        println!("no sequences found.");
        exit(1);
    }

    println!("seq_id\tpos\twindow\tTm_basic_C\tTm_salt_C");

    for (id, seq_raw) in entries {
        let seq: String = seq_raw
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| c.to_ascii_uppercase())
            .collect();

        if seq.len() < W {
            eprintln!("{} length < {}. skipping.", id, W);
            continue;
        }

        if let Some(bad) = seq.chars().find(|&c| !matches!(c, 'A' | 'C' | 'G' | 'T')) {
            println!("invalid char '{}' in seq {}", bad, id);
            exit(1);
        }

        let bytes = seq.as_bytes();
        let mut cnt = [0usize; 4];

        for i in 0..W {
            cnt[idx(bytes[i])] += 1;
        }

        for i in 0..=bytes.len() - W {
            let a = cnt[0] as f64;
            let c = cnt[1] as f64;
            let g = cnt[2] as f64;
            let t = cnt[3] as f64;
            let len = W as f64;

            let tm_basic = 4.0 * (g + c) + 2.0 * (a + t);
            let gc_percent = 100.0 * (g + c) / len;
            let tm_salt = 81.5 + 16.6 * NA_MOLAR.log10() + 0.41 * gc_percent - 600.0 / len;

            let window = &seq[i..i + W];
            println!("{}\t{}\t{}\t{:.2}\t{:.2}", id, i + 1, window, tm_basic, tm_salt);

            if i + W < bytes.len() {
                let out = bytes[i];
                let inn = bytes[i + W];
                cnt[idx(out)] -= 1;
                cnt[idx(inn)] += 1;
            }
        }
    }
}

fn idx(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => unreachable!(),
    }
}

fn parse_fasta(s: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    let mut id: Option<String> = None;
    let mut seq = String::new();

    for line in s.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if let Some(cur_id) = id.take() {
                out.push((cur_id, seq.clone()));
                seq.clear();
            }
            id = Some(line[1..].trim().to_string());
        } else {
            seq.push_str(line);
        }
    }
    if let Some(cur_id) = id {
        out.push((cur_id, seq));
    } else if !seq.is_empty() {
        out.push(("seq1".to_string(), seq));
    }
    out
}
