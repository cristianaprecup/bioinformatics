use std::env;
use std::fs;
use std::io::{self, Read};

#[derive(Debug)]
struct TandemRepeat {
    start: usize,
    k: usize,
    motif: String,
    repeats: usize,
}

fn extract_sequence(raw: &str) -> String {
    let mut seq = String::new();
    for line in raw.lines() {
        if line.starts_with('>') {
            continue;
        }
        seq.push_str(line);
    }
    let mut cleaned = String::with_capacity(seq.len());
    for ch in seq.chars() {
        match ch.to_ascii_uppercase() {
            'A' | 'C' | 'G' | 'T' | 'N' => cleaned.push(ch.to_ascii_uppercase()),
            _ => {}
        }
    }
    cleaned
}

fn read_sequence_from_path(path: &str) -> io::Result<String> {
    Ok(extract_sequence(&fs::read_to_string(path)?))
}

fn read_sequence_from_stdin() -> io::Result<String> {
    let mut raw = String::new();
    io::stdin().read_to_string(&mut raw)?;
    Ok(extract_sequence(&raw))
}

fn find_tandem_repeats(seq: &str, k_min: usize, k_max: usize) -> Vec<TandemRepeat> {
    let b = seq.as_bytes();
    let n = b.len();
    let mut hits: Vec<TandemRepeat> = Vec::new();

    let mut i = 0usize;
    while i + k_min * 2 <= n {
        let mut advanced = false;

        for k in k_min..=k_max {
            if i + 2 * k > n {
                break;
            }

            let motif = &b[i..i + k];

            let mut r = 1usize;
            while i + k * (r + 1) <= n && &b[i + k * r..i + k * (r + 1)] == motif {
                r += 1;
            }

            if r >= 2 {
                hits.push(TandemRepeat {
                    start: i,
                    k,
                    motif: String::from_utf8(motif.to_vec()).unwrap(),
                    repeats: r,
                });
                i += k * r;
                advanced = true;
                break;
            }
        }

        if !advanced {
            i += 1;
        }
    }

    hits
}
use plotters::prelude::*;

fn plot_histogram_png(
    filename: &str,
    title: &str,
    x_labels: &[String],
    values: &[u32],
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    let n = x_labels.len();
    let root = BitMapBackend::new(filename, (1000, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    if n == 0 {
        let mut chart = ChartBuilder::on(&root)
            .caption(title, ("sans-serif", 28))
            .margin(20)
            .build_cartesian_2d(0..1, 0u32..1)?;
        chart.configure_mesh().draw()?;
        root.present()?;
        return Ok(());
    }

    let ymax = values.iter().copied().max().unwrap_or(0).max(1);
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(70)
        .y_label_area_size(60)
        .build_cartesian_2d(0..(n as i32), 0u32..(ymax + (ymax / 5).max(1)))?;

    chart
        .configure_mesh()
        .x_labels(n)
        .x_label_formatter(&|x| {
            let i = (*x as usize).min(n.saturating_sub(1));
            x_labels[i].clone()
        })
        .y_desc("frequency")
        .x_desc("Category")
        .label_style(("sans-serif", 14))
        .axis_desc_style(("sans-serif", 18))
        .draw()?;

    chart.draw_series((0..n).map(|i| {
        let v = values[i];
        let x0 = i as i32;
        let x1 = x0 + 1;
        Rectangle::new([(x0, 0u32), (x1, v)], BLUE.filled())
    }))?;

    chart.draw_series(
        (0..n)
            .filter(|&i| values[i] > 0)
            .map(|i| {
                let v = values[i];
                let xm = i as i32; 
                Text::new(v.to_string(), (xm, v + (ymax / 50).max(1)), ("sans-serif", 14))
            }),
    )?;

    root.present()?;
    Ok(())
}


fn plot_frequencies(
    hits: &[TandemRepeat],
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::BTreeMap;

    let mut by_k: BTreeMap<usize, u32> = BTreeMap::new();
    for h in hits {
        *by_k.entry(h.k).or_insert(0) += 1;
    }
    let xk: Vec<String> = by_k.keys().map(|k| k.to_string()).collect();
    let yk: Vec<u32> = by_k.values().copied().collect();
    plot_histogram_png(
        "repeats_by_k.png",
        "Number of repeats per motif length k",
        &xk,
        &yk,
    )?;

    let mut by_r: BTreeMap<usize, u32> = BTreeMap::new();
    for h in hits {
        *by_r.entry(h.repeats).or_insert(0) += 1;
    }
    let xr: Vec<String> = by_r.keys().map(|r| r.to_string()).collect();
    let yr: Vec<u32> = by_r.values().copied().collect();
    plot_histogram_png(
        "repeats_by_repeats.png",
        "Number of tandem repeats per multiplicity r",
        &xr,
        &yr,
    )?;
    
    Ok(())
}

fn plot_motif_barchart_r3(
    out_basename: &str,   
    display_name: &str,   
    hits: &[TandemRepeat],
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::HashMap;
    use plotters::prelude::*;
    // Numără DOAR motivele cu r = 3
    let mut count_by_motif: HashMap<String, u32> = HashMap::new();
    for h in hits.iter().filter(|h| h.repeats == 3) {
        *count_by_motif.entry(h.motif.clone()).or_insert(0) += 1;
    }
    if count_by_motif.is_empty() {
        return Ok(());
    }

    let mut items: Vec<(String, u32)> = count_by_motif.into_iter().collect();
    items.sort_by_key(|(_, c)| std::cmp::Reverse(*c));

    let labels: Vec<String> = items.iter().map(|(m, _)| m.clone()).collect();
    let values: Vec<u32>   = items.iter().map(|(_, c)| *c).collect();
    let xmax = values.iter().copied().max().unwrap_or(1);

    let h = 120 + (labels.len() as u32) * 24;
    let w = 1100u32;
    let outfile = format!("{}_r3_motifs.png", out_basename);
    let root = BitMapBackend::new(&outfile, (w, h)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(format!("Tandem repeats r = 3 in {}", display_name), ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(50)
        .y_label_area_size(240)
        .build_cartesian_2d(0u32..(xmax + (xmax / 5).max(1)), 0..(labels.len() as i32))?;

    chart
        .configure_mesh()
        .y_labels(labels.len())
        .y_label_formatter(&|i| {
            let idx = (*i as usize).min(labels.len().saturating_sub(1));
            labels[idx].clone()
        })
        .x_desc("Frequency")
        .y_desc("Motif (3–10 bases, r = 3)")
        .label_style(("sans-serif", 14))
        .axis_desc_style(("sans-serif", 18))
        .draw()?;

    chart.draw_series((0..labels.len()).map(|i| {
        let v = values[i];
        let y0 = i as i32;
        let y1 = y0 + 1;
        Rectangle::new([(0u32, y0), (v, y1)], BLUE.filled())
    }))?;

    chart.draw_series((0..labels.len()).map(|i| {
        let v = values[i];
        let y = i as i32;
        Text::new(v.to_string(), (v, y), ("sans-serif", 14))
    }))?;

    root.present()?;
    Ok(())
}

fn file_stem_from_arg(arg: &str) -> String {
    std::path::Path::new(arg)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("seq")
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}


fn main() -> std::io::Result<()> {
    use std::env;
    let args: Vec<String> = env::args().collect();

    const K_MIN: usize = 3;
    const K_MAX: usize = 10;

    let seq = if args.len() > 1 {
        read_sequence_from_path(&args[1])?
    } else {
        read_sequence_from_stdin()?
    };

    
    let hits = find_tandem_repeats(&seq, K_MIN, K_MAX);
    
    let basename = if args.len() > 1 {
    file_stem_from_arg(&args[1])
} else {
    "seq".to_string()
};

if let Err(e) = plot_motif_barchart_r3(&basename, &basename, &hits) {
    eprintln!("could not plot r=3 motif chart for {}: {}", basename, e);
}


    println!("start ,k, seq, repeats, total_len");
    for tr in find_tandem_repeats(&seq, K_MIN, K_MAX) {
        let total_len = tr.k * tr.repeats;
        println!(
            "start: {}, k: {}, sequence: {}, repeats: {}, total_len: {}",
            tr.start, tr.k, tr.motif, tr.repeats, total_len
        );
    }

    if let Err(e) = plot_frequencies(&hits) {
        println!("could not generate plots, {e}");
    }

    Ok(())
}
