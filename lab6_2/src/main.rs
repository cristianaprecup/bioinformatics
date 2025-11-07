use std::fs;
use std::path::Path;
use rand::{rngs::StdRng, SeedableRng};
use std::cmp::{min};

#[derive(Debug)]
struct Fragment {
    start: usize,
    end: usize,
    length: usize,
}

fn read_fasta_sequence<P: AsRef<Path>>(path: P) -> Result<String, String> {
    let content = fs::read_to_string(&path).map_err(|e| format!("cant read: {}", e))?;
    let mut seq = String::new();
    for line in content.lines() {
        let line = line.trim();
        if !line.starts_with('>') && !line.is_empty() {
            seq.push_str(line);
        }
    }
    Ok(seq.to_uppercase())
}

fn digest_ecori(sequence: &str) -> Vec<Fragment> {
    let pattern = "GAATTC";
    let mut fragments = Vec::new();
    let mut start = 0;
    let mut search_pos = 0;

    while let Some(pos) = sequence[search_pos..].find(pattern) {
        let cut_site = search_pos + pos + 1; // G^AATTC
        fragments.push(Fragment {
            start,
            end: cut_site,
            length: cut_site - start,
        });
        start = cut_site;
        search_pos = cut_site;
    }

    if start < sequence.len() {
        fragments.push(Fragment {
            start,
            end: sequence.len(),
            length: sequence.len() - start,
        });
    }

    fragments
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

fn draw_gel(fragments: &[Fragment], title: &str) -> String {
    let gel_height = 34usize;
    let gel_width = 30usize;
    let mut canvas = vec![vec![' '; gel_width]; gel_height];

    for y in 0..gel_height {
        canvas[y][0] = '|';
        canvas[y][gel_width - 1] = '|';
    }
    for x in 0..gel_width {
        canvas[0][x] = if x == 0 || x == gel_width - 1 { '+' } else { '-' };
    }

    let min_bp = fragments.iter().map(|f| f.length).min().unwrap_or(100);
    let max_bp = fragments.iter().map(|f| f.length).max().unwrap_or(3000);

    for frag in fragments {
        let y = min(
            gel_height - 1,
            len_to_position(frag.length, min_bp, max_bp, gel_height),
        );
        for x in 2..gel_width - 2 {
            canvas[y][x] = '=';
        }
    }

    let mut out = format!("Gel electrophoresis: {}\n\n", title);
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

    let paths = fs::read_dir(".").expect("cant read current dir");
    let fasta_files: Vec<_> = paths
        .filter_map(|p| {
            let path = p.ok()?.path();
            if path.extension()? == "fna" {
                Some(path)
            } else {
                None
            }
        })
        .collect();


    for file in fasta_files {
        let name = file.file_name().unwrap().to_string_lossy().to_string();
        println!("\nAnalyze: {}", name);

        let seq = read_fasta_sequence(&file).expect("cant read sequence");
        let fragments = digest_ecori(&seq);

        println!("Number of fragments: {}", fragments.len());
        println!(
            "Fragment lengths: {:?}",
            fragments.iter().map(|f| f.length).collect::<Vec<_>>()
        );

        let gel = draw_gel(&fragments, &name);
        println!("{}", gel);
    }
}
