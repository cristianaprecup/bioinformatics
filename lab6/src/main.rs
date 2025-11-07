use rand::{Rng, SeedableRng, rngs::StdRng};
use std::cmp::{max, min};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
struct Fragment {
    start: usize,
    len: usize,
    seq: String,
}

fn read_fasta_first_record<P: AsRef<Path>>(path: P) -> Result<(String, String), String> {
    let content = fs::read_to_string(&path).map_err(|e| format!("cant read file: {}", e))?;
    let mut header = String::new();
    let mut seq = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if seq.is_empty() {
                header = line[1..].to_string();
            } else {
                break;
            }
        } else if !line.starts_with(';') {
            for ch in line.chars() {
                let up = ch.to_ascii_uppercase();
                if matches!(up, 'A' | 'C' | 'G' | 'T' | 'N') {
                    seq.push(up);
                }
            }
        }
    }

    Ok((header, seq))
}

fn len_to_position(len: usize, min_bp: usize, max_bp: usize, gel_height: usize) -> usize {
    let eps = 1e-9_f64;
    let lmin = (min_bp as f64).log10();
    let lmax = (max_bp as f64).log10();
    let lcur = (len as f64).log10();
    let norm = (lmax - lcur) / ((lmax - lmin) + eps);
    let top_margin = 2usize;
    let bottom_margin = 1usize;
    let usable = gel_height.saturating_sub(top_margin + bottom_margin);
    top_margin + (norm * usable as f64).round() as usize
}

fn draw_gel(fragments: &[Fragment]) -> String {
    let gel_height = 34usize;
    let gel_width = 30usize;
    let mut canvas = vec![vec![' '; gel_width]; gel_height];

    for y in 0..gel_height {
        canvas[y][0] = '|';
        canvas[y][gel_width - 1] = '|';
    }
    for x in 0..gel_width {
        canvas[0][x] = if x == 0 || x == gel_width - 1 {
            '+'
        } else {
            '-'
        };
    }

    let min_bp = fragments.iter().map(|f| f.len).min().unwrap_or(100);
    let max_bp = fragments.iter().map(|f| f.len).max().unwrap_or(3000);

    for frag in fragments {
        let y = min(
            gel_height - 1,
            len_to_position(frag.len, min_bp, max_bp, gel_height),
        );
        for x in 2..gel_width - 2 {
            canvas[y][x] = '=';
        }
    }

    let mut out = String::new();
    for row in canvas {
        for ch in row {
            out.push(ch);
        }
        out.push('\n');
    }
    out
}

fn main() {
    let mut rng = StdRng::seed_from_u64(42);

    let path = std::env::args().nth(1).expect("lab6 <path_to_fasta>");

    let (header, genome) = read_fasta_first_record(path).expect("error reading FASTA");

    let mut fragments = Vec::new();
    for _ in 0..10 {
        let max_len_here = max(100, min(3000, genome.len()));
        let len = rng.gen_range(100..=max_len_here);
        let start_max = genome.len().saturating_sub(len);
        let start = if start_max == 0 {
            0
        } else {
            rng.gen_range(0..=start_max)
        };
        let seq = genome[start..start + len].to_string();
        fragments.push(Fragment { start, len, seq });
    }

    println!("Header FASTA: {}", header);
    println!("FASTA length: {} nt", genome.len());
    println!("Extracted fragments: {}", fragments.len());
    println!();

    println!("{:<6} {:<8} {:<8}", "Idx", "Start", "Length(bp)");
    for (i, f) in fragments.iter().enumerate() {
        println!("{:<6} {:<8} {:<8}", i, f.start, f.len);
    }

    println!("\nGel electrophoresis (ASCII):\n");
    let gel = draw_gel(&fragments);
    println!("{}", gel);

    for (i, f) in fragments.iter().take(10).enumerate() {
        let preview: String = f.seq.chars().take(40).collect();
        println!(
            "Frag {} preview: {}{}",
            i,
            preview,
            if f.len > 40 { "..." } else { "" }
        );
    }
}
