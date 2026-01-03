use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let s = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG";
    let window_size = 30;

    println!("Sequence: {}", s);
    println!("Total Length: {}", s.len());
    println!("Window Size: {}", window_size);
    println!("---");

    let mut cg_values = Vec::new();
    let mut kappa_values = Vec::new();
    let mut points = Vec::new();

    let windows: Vec<&str> = s
        .as_bytes()
        .windows(window_size)
        .map(|w| std::str::from_utf8(w).unwrap())
        .collect();

    for (i, window) in windows.iter().enumerate() {
        let cg = calc_cg_content(window);
        
        let kappa = calc_kappa_ic(window);

        cg_values.push(cg);
        kappa_values.push(kappa);
        points.push((cg, kappa));

        if i == 0 {
             println!("Window 0: {}", window);
             println!("CG%: {:.2}", cg);
             println!("Kappa IC: {:.2}", kappa);
        }
    }

    let avg_cg: f64 = cg_values.iter().sum::<f64>() / cg_values.len() as f64;
    let avg_kappa: f64 = kappa_values.iter().sum::<f64>() / kappa_values.len() as f64;

    println!("---");
    println!("Calculated Average CG%:     {:.2}", avg_cg);
    println!("Calculated Average Kappa IC: {:.2}", avg_kappa);

    plot_pattern(&points, "pattern_chart.png")?;
    
    plot_center(avg_cg, avg_kappa, "center_chart.png")?;

    Ok(())
}

fn calc_cg_content(window: &str) -> f64 {
    let mut count = 0;
    for c in window.chars() {
        if c == 'C' || c == 'G' {
            count += 1;
        }
    }
    (count as f64 / window.len() as f64) * 100.0
}

fn calc_kappa_ic(window: &str) -> f64 {
    let n = window.len();
    let chars: Vec<char> = window.chars().collect();
    let mut total_probability = 0.0;

    let max_u = n - 1; 

    for u in 1..=max_u {
        let mut matches = 0;
        let suffix_len = n - u;
        
        for i in 0..suffix_len {
            if chars[i] == chars[i + u] {
                matches += 1;
            }
        }

        if suffix_len > 0 {
            total_probability += (matches as f64 / suffix_len as f64) * 100.0;
        }
    }

    total_probability / max_u as f64
}

fn plot_pattern(points: &[(f64, f64)], filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("DNA Pattern (C+G% vs Kappa IC)", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..100f64, 0f64..100f64)?;

    chart.configure_mesh()
        .x_desc("C+G %")
        .y_desc("Kappa IC")
        .draw()?;

    chart.draw_series(
        points.iter().map(|(x, y)| Circle::new((*x, *y), 3, BLUE.filled()))
    )?;

    Ok(())
}

fn plot_center(avg_cg: f64, avg_kappa: f64, filename: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(filename, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Center of Weight", ("sans-serif", 30).into_font())
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..100f64, 0f64..100f64)?;

    chart.configure_mesh()
        .x_desc("Average C+G %")
        .y_desc("Average Kappa IC")
        .draw()?;

    chart.draw_series(std::iter::once(
        EmptyElement::at((avg_cg, avg_kappa))
            + Circle::new((0, 0), 6, RED.filled())
            + Text::new(
                format!("({:.2}, {:.2})", avg_cg, avg_kappa),
                (10, 0),
                ("sans-serif", 15).into_font(),
            ),
    ))?;

    Ok(())
}