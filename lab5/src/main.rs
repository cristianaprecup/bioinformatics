use rand::{Rng, SeedableRng, rngs::StdRng};
use std::collections::{HashMap, HashSet};
use std::env;
use std::error::Error;
use std::fs;

fn arg_k() -> Option<usize> {
    let mut it = env::args().skip(2);
    while let Some(tok) = it.next() {
        if tok == "--k" {
            if let Some(v) = it.next() {
                if let Ok(k) = v.parse() {
                    return Some(k);
                }
            }
        } else if let Some(v) = tok.strip_prefix("--k=") {
            if let Ok(k) = v.parse() {
                return Some(k);
            }
        }
    }
    None
}

fn arg_seed() -> Option<u64> {
    let mut it = env::args().skip(2);
    while let Some(tok) = it.next() {
        if tok == "--seed" {
            if let Some(v) = it.next() {
                if let Ok(s) = v.parse() {
                    return Some(s);
                }
            }
        } else if let Some(v) = tok.strip_prefix("--seed=") {
            if let Ok(s) = v.parse() {
                return Some(s);
            }
        }
    }
    None
}

fn arg_record() -> Option<String> {
    let mut it = env::args().skip(2);
    while let Some(tok) = it.next() {
        if tok == "--record" {
            if let Some(v) = it.next() {
                return Some(v);
            }
        } else if let Some(v) = tok.strip_prefix("--record=") {
            return Some(v.to_string());
        }
    }
    None
}

fn read_fasta_records(path: &str) -> Result<Vec<(String, String)>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let mut records: Vec<(String, String)> = Vec::new();
    let mut header = String::new();
    let mut seq = String::new();

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if !header.is_empty() || !seq.is_empty() {
                if !seq.is_empty() {
                    records.push((header.clone(), seq.clone()));
                }
                seq.clear();
            }
            header = line[1..].to_string();
        } else {
            for ch in line.chars() {
                let up = match ch {
                    'a' | 'A' => 'A',
                    'c' | 'C' => 'C',
                    'g' | 'G' => 'G',
                    't' | 'T' | 'u' | 'U' => 'T',
                    'n' | 'N' => 'N',
                    _ => continue,
                };
                seq.push(up);
            }
        }
    }
    if !header.is_empty() && !seq.is_empty() {
        records.push((header, seq));
    }

    if records.is_empty() {
        return Err("No FASTA records found".into());
    }
    Ok(records)
}

fn sample_reads(seq: &str, n: usize, len_min: usize, len_max: usize, seed: u64) -> Vec<String> {
    assert!(len_min > 0 && len_min <= len_max);
    assert!(seq.len() >= len_max, "Sequence must be longer than reads");

    let mut rng = StdRng::seed_from_u64(seed);
    let mut reads = Vec::with_capacity(n);
    for _ in 0..n {
        let len = rng.random_range(len_min..=len_max);
        let start = rng.random_range(0..=seq.len() - len);
        reads.push(seq[start..start + len].to_string());
    }
    reads
}

fn assemble_debruijn(reads: &[String], k: usize) -> String {
    assert!(k >= 2, "k trebuie să fie ≥ 2");
    let min_read = reads.iter().map(|r| r.len()).min().expect("fără reads");
    assert!(k <= min_read, "k trebuie ≤ lungimea minimă a read-ului");

    let mut node_id: HashMap<String, usize> = HashMap::new();
    let mut id_node: Vec<String> = Vec::new();
    let mut next_id = 0usize;

    let mut adj: Vec<Vec<(usize, u8)>> = Vec::new();
    let mut indeg: Vec<usize> = Vec::new();
    let mut outdeg: Vec<usize> = Vec::new();

    let intern = |s: &str,
                  node_id: &mut HashMap<String, usize>,
                  id_node: &mut Vec<String>,
                  adj: &mut Vec<Vec<(usize, u8)>>,
                  indeg: &mut Vec<usize>,
                  outdeg: &mut Vec<usize>,
                  next_id: &mut usize| {
        if let Some(&id) = node_id.get(s) {
            id
        } else {
            let id = *next_id;
            *next_id += 1;
            node_id.insert(s.to_string(), id);
            id_node.push(s.to_string());
            adj.push(Vec::new());
            indeg.push(0);
            outdeg.push(0);
            id
        }
    };

    let mut seen_kmers: HashSet<String> = HashSet::new();

    for read in reads {
        if read.len() < k {
            continue;
        }
        for i in 0..=read.len() - k {
            let kmer = &read[i..i + k];
            if !seen_kmers.insert(kmer.to_string()) {
                continue;
            }
            let prefix = &kmer[..k - 1];
            let suffix = &kmer[1..];
            let p = intern(
                prefix,
                &mut node_id,
                &mut id_node,
                &mut adj,
                &mut indeg,
                &mut outdeg,
                &mut next_id,
            );
            let s = intern(
                suffix,
                &mut node_id,
                &mut id_node,
                &mut adj,
                &mut indeg,
                &mut outdeg,
                &mut next_id,
            );
            let ch = kmer.as_bytes()[k - 1];
            adj[p].push((s, ch));
            outdeg[p] += 1;
            indeg[s] += 1;
        }
    }

    if next_id == 0 {
        return String::new();
    }

    let mut start = 0usize;
    for v in 0..next_id {
        if outdeg[v] == indeg[v] + 1 {
            start = v;
            break;
        }
        if outdeg[start] == 0 && indeg[start] == 0 && (outdeg[v] + indeg[v] > 0) {
            start = v;
        }
    }

    let mut stack: Vec<usize> = vec![start];
    let mut edge_char_stack: Vec<u8> = Vec::new();
    let mut out_chars: Vec<u8> = Vec::new();

    while let Some(&v) = stack.last() {
        if let Some((next, ch)) = adj[v].pop() {
            stack.push(next);
            edge_char_stack.push(ch);
        } else {
            stack.pop();
            if let Some(ch) = edge_char_stack.pop() {
                out_chars.push(ch);
            }
        }
    }
    out_chars.reverse();

    if out_chars.is_empty() {
        return String::new();
    }

    let mut seq = id_node[start].clone().into_bytes();
    seq.extend(out_chars);
    String::from_utf8(seq).unwrap()
}

fn choose_k(reads: &[String]) -> usize {
    let min_read = reads.iter().map(|r| r.len()).min().unwrap_or(100);
    let mut k = 51usize;
    if k >= min_read {
        k = std::cmp::max(21, min_read.saturating_sub(1));
    }
    if k % 2 == 0 {
        k -= 1;
    }
    k
}

fn revcomp(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars().rev() {
        out.push(match ch {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        });
    }
    out
}

fn kmers(s: &str, k: usize) -> HashSet<String> {
    let mut set = HashSet::new();
    if s.len() >= k {
        for i in 0..=s.len() - k {
            set.insert(s[i..i + k].to_string());
        }
    }
    set
}

fn fraction_kmers_in(a: &str, b: &str, k: usize) -> f64 {
    if a.len() < k {
        return 0.0;
    }
    let aa = kmers(a, k);
    let bb = kmers(b, k);
    if aa.is_empty() {
        return 0.0;
    }
    let hit = aa.iter().filter(|x| bb.contains(*x)).count();
    hit as f64 / aa.len() as f64
}

fn print_wrapped(title: &str, seq: &str, width: usize) {
    println!("{}", title);
    if seq.is_empty() {
        println!("[gol]\n");
        return;
    }
    for chunk in seq.as_bytes().chunks(width) {
        println!("{}", std::str::from_utf8(chunk).unwrap());
    }
    println!();
}

fn main() -> Result<(), Box<dyn Error>> {
    let path = env::args()
        .nth(1)
        .expect("Provide the path to the .fasta file as arg");

    let seed = match arg_seed() {
        Some(s) => {
            println!("Chosen seed: {}", s);
            s
        }
        None => {
            let s: u64 = rand::thread_rng().random();
            println!("chosen seed: {}", s);
            s
        }
    };
    let mut rng_for_choose = StdRng::seed_from_u64(seed);

    let records = read_fasta_records(&path)?;
    println!("FASTAs found: {}", records.len());

    let original_idx = match arg_record() {
        Some(s) if s.eq_ignore_ascii_case("random") => {
            rng_for_choose.random_range(0..records.len())
        }
        Some(s) => s
            .parse::<usize>()
            .ok()
            .map(|i| if i < records.len() { i } else { 0 })
            .unwrap_or(0),
        None => 0,
    };

    let (header, original) = &records[original_idx];
    println!("Record chosen: #{} | {}", original_idx, header);
    println!("Read {} bases from {}", original.len(), path);

    let seed_reads: u64 = rng_for_choose.random();
    let reads = sample_reads(original.as_str(), 2000, 100, 150, seed_reads);
    let avg_len: f64 = reads.iter().map(|r| r.len()).sum::<usize>() as f64 / reads.len() as f64;
    let cov = (reads.len() as f64) * avg_len / (original.len().max(1) as f64);
    println!(
        "Generated {} reads. Average length ≈ {:.1}. Estimated coverage ≈ {:.1}",
        reads.len(),
        avg_len,
        cov
    );

    let k = arg_k().unwrap_or_else(|| choose_k(&reads));
    let assembled = assemble_debruijn(&reads, k);

    println!("k used: {}", k);
    println!(
        "Original length: {} | Assembled length: {}",
        original.len(),
        assembled.len()
    );
    println!();

    print_wrapped(">>> ORIGINAL", original.as_str(), 80);
    print_wrapped(">>> RECONSTRUCTION", assembled.as_str(), 80);

    if assembled.is_empty() {
        println!("Empty assembly.");
        return Ok(());
    }

    if let Some(pos) = original.as_str().find(assembled.as_str()) {
        println!(
            "Original appears as a substring in the contig at position {}.",
            pos
        );
        if assembled.len() == original.len() {
            println!("Same reconstruction.");
        } else {
            println!("Partial reconstruction.");
        }
        return Ok(());
    }

    let assembled_rc = revcomp(assembled.as_str());
    let original_rc = revcomp(original.as_str());

    let frac_fwd = fraction_kmers_in(assembled.as_str(), original.as_str(), k);
    let frac_rc = fraction_kmers_in(assembled_rc.as_str(), original.as_str(), k);
    let frac = frac_fwd.max(frac_rc);

    println!(
        "With k = {}, fraction of k-mers from the assembled contig found in the original (or its reverse complement) = {:.2}%.",
        100.0 * frac,
        k
    );

    Ok(())
}
