use plotters::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs;
use std::path::Path;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        std::process::exit(1);
    }

    let covid_seq = read_fasta(&args[1]);
    let flu_seq = read_fasta(&args[2]);

    let covid_counts = codon_counts(&covid_seq);
    let flu_counts = codon_counts(&flu_seq);

    let covid_freq = counts_to_freq(&covid_counts);
    let flu_freq = counts_to_freq(&flu_counts);

    plot_top10(
        &covid_freq,
        "Top 10 codons COVID-19 frequency",
        "covid_top10.png",
    )?;

    plot_top10(
        &flu_freq,
        "Top 10 codons Influenza frequency",
        "influenza_top10.png",
    )?;

    plot_compare_top(
        &covid_freq,
        &flu_freq,
        "Comparison of COVID-19 vs Influenza codon frequencies",
        "compare_top10.png",
    )?;

    println!();
    println!("Top 3 aminoacis COVID-19:");
    print_top3_aa(&covid_counts);

    println!();
    println!("Top 3 aminoacids Influenza:");
    print_top3_aa(&flu_counts);

    Ok(())
}

fn read_fasta<P: AsRef<Path>>(path: P) -> String {
    let raw = fs::read_to_string(path).expect("Cannot read file");
    let seq: String = raw
        .lines()
        .filter(|l| !l.starts_with('>'))
        .collect::<Vec<_>>()
        .join("");

    seq.to_uppercase()
        .replace('T', "U")
        .chars()
        .filter(|c| matches!(c, 'A' | 'C' | 'G' | 'U'))
        .collect()
}

fn codon_counts(rna: &str) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    let bytes = rna.as_bytes();
    let mut i = 0;
    while i + 3 <= bytes.len() {
        let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap().to_string();
        *map.entry(codon).or_insert(0) += 1;
        i += 3;
    }
    map
}

fn counts_to_freq(counts: &HashMap<String, usize>) -> HashMap<String, f64> {
    let total: usize = counts.values().sum();
    if total == 0 {
        return HashMap::new();
    }
    counts
        .iter()
        .map(|(k, &v)| (k.clone(), v as f64 / total as f64))
        .collect()
}

fn top_n(map: &HashMap<String, f64>, n: usize) -> Vec<(String, f64)> {
    let mut v: Vec<(String, f64)> = map.iter().map(|(k, &v)| (k.clone(), v)).collect();
    v.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    v.truncate(n);
    v
}
fn plot_top10(freqs: &HashMap<String, f64>, title: &str, out: &str) -> Result<(), Box<dyn Error>> {
    let data = top_n(freqs, 10);
    let labels: Vec<String> = data.iter().map(|(c, _)| c.clone()).collect();
    let values: Vec<f64> = data.iter().map(|(_, f)| f * 100.0).collect();

    if labels.is_empty() {
        eprintln!("did not find any codons to plot for {}", title);
        return Ok(());
    }

    let root = BitMapBackend::new(out, (1000, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let ymax = values.iter().cloned().fold(0.0_f64, f64::max).max(1.0) * 1.15;
    let n = labels.len();

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(0..n as i32, 0.0..ymax)?;

    chart
        .configure_mesh()
        .x_labels(n)
        .x_label_formatter(&|x| {
            let i = (*x as usize).min(n.saturating_sub(1));
            labels[i].clone()
        })
        .y_label_formatter(&|y| format!("{y:.1}%"))
        .axis_desc_style(("sans-serif", 16))
        .draw()?;

    chart.draw_series(values.iter().enumerate().map(|(i, v)| {
        let x0 = i as i32;
        let x1 = x0 + 1;
        let mut bar = Rectangle::new([(x0, 0.0), (x1, *v)], ShapeStyle::from(&BLUE).filled());
        bar.set_margin(5, 5, 0, 0);
        bar
    }))?;

    Ok(())
}

fn plot_compare_top(
    a: &HashMap<String, f64>,
    b: &HashMap<String, f64>,
    title: &str,
    out: &str,
) -> Result<(), Box<dyn Error>> {
    let top_a: HashSet<String> = top_n(a, 10).into_iter().map(|x| x.0).collect();
    let top_b: HashSet<String> = top_n(b, 10).into_iter().map(|x| x.0).collect();

    let union: HashSet<String> = top_a.union(&top_b).cloned().collect();

    let mut combined: Vec<(String, f64)> = union
        .into_iter()
        .map(|c| {
            let fa = *a.get(&c).unwrap_or(&0.0);
            let fb = *b.get(&c).unwrap_or(&0.0);
            (c, (fa + fb) / 2.0)
        })
        .collect();

    combined.sort_by(|x, y| y.1.partial_cmp(&x.1).unwrap());
    combined.truncate(10);

    let labels: Vec<String> = combined.iter().map(|x| x.0.clone()).collect();

    let vals_a: Vec<f64> = labels
        .iter()
        .map(|c| a.get(c).copied().unwrap_or(0.0) * 100.0)
        .collect();
    let vals_b: Vec<f64> = labels
        .iter()
        .map(|c| b.get(c).copied().unwrap_or(0.0) * 100.0)
        .collect();

    let root = BitMapBackend::new(out, (1200, 650)).into_drawing_area();
    root.fill(&WHITE)?;

    let ymax = vals_a
        .iter()
        .chain(vals_b.iter())
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(1.0)
        * 1.20;

    let n = labels.len();
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(60)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..n as f64, 0.0..ymax)?;

    chart
        .configure_mesh()
        .x_labels(n)
        .x_label_formatter(&|x| {
            let i = x.floor().clamp(0.0, (n.saturating_sub(1)) as f64) as usize;
            labels[i].clone()
        })
        .y_label_formatter(&|y| format!("{y:.1}%"))
        .axis_desc_style(("sans-serif", 16))
        .draw()?;

    for i in 0..n {
        let xa0 = i as f64 + 0.08;
        let xa1 = i as f64 + 0.45;
        let xb0 = i as f64 + 0.55;
        let xb1 = i as f64 + 0.92;

        chart.draw_series(std::iter::once(Rectangle::new(
            [(xa0, 0.0), (xa1, vals_a[i])],
            ShapeStyle::from(&BLUE).filled(),
        )))?;

        chart.draw_series(std::iter::once(Rectangle::new(
            [(xb0, 0.0), (xb1, vals_b[i])],
            ShapeStyle::from(&RED).filled(),
        )))?;
    }

    chart
        .draw_series(std::iter::empty::<Rectangle<_>>())?
        .label("COVID-19")
        .legend(|(x, y)| Rectangle::new([(x, y - 5), (x + 20, y + 5)], BLUE.filled()));

    chart
        .draw_series(std::iter::empty::<Rectangle<_>>())?
        .label("Influenza")
        .legend(|(x, y)| Rectangle::new([(x, y - 5), (x + 20, y + 5)], RED.filled()));

    chart
        .configure_series_labels()
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

fn print_top3_aa(counts: &HashMap<String, usize>) {
    let mut aa_counts: HashMap<&'static str, usize> = HashMap::new();
    for (codon, &cnt) in counts {
        match codon_to_aa(codon) {
            Some(name) if name != "Stop" => {
                *aa_counts.entry(name).or_insert(0) += cnt;
            }
            _ => {}
        }
    }
    let mut v: Vec<(&str, usize)> = aa_counts.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1));
    v.truncate(3);

    for (i, (aa, n)) in v.iter().enumerate() {
        println!("{}. {}  {}", i + 1, aa, n);
    }
}

fn codon_to_aa(c: &str) -> Option<&'static str> {
    Some(match c {
        "UUU" | "UUC" => "Phe",
        "UUA" | "UUG" | "CUU" | "CUC" | "CUA" | "CUG" => "Leu",
        "AUU" | "AUC" | "AUA" => "Ile",
        "AUG" => "Met",
        "GUU" | "GUC" | "GUA" | "GUG" => "Val",
        "UCU" | "UCC" | "UCA" | "UCG" | "AGU" | "AGC" => "Ser",
        "CCU" | "CCC" | "CCA" | "CCG" => "Pro",
        "ACU" | "ACC" | "ACA" | "ACG" => "Thr",
        "GCU" | "GCC" | "GCA" | "GCG" => "Ala",
        "UAU" | "UAC" => "Tyr",
        "CAU" | "CAC" => "His",
        "CAA" | "CAG" => "Gln",
        "AAU" | "AAC" => "Asn",
        "AAA" | "AAG" => "Lys",
        "GAU" | "GAC" => "Asp",
        "GAA" | "GAG" => "Glu",
        "UGU" | "UGC" => "Cys",
        "UGG" => "Trp",
        "CGU" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => "Arg",
        "GGU" | "GGC" | "GGA" | "GGG" => "Gly",
        "UAA" | "UAG" | "UGA" => "Stop",
        _ => return None,
    })
}
