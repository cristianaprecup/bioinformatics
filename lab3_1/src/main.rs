use std::io;

fn main() {
    let mut seq = String::new();
    io::stdin().read_line(&mut seq).expect("error at reading stdin");

    let seq = seq.trim().to_uppercase();
    if seq.is_empty() {
        println!("error: empty seq.");
        std::process::exit(1);
    }

    let mut a = 0usize;
    let mut c = 0usize;
    let mut g = 0usize;
    let mut t = 0usize;

    for ch in seq.chars() {
        match ch {
            'A' => a += 1,
            'C' => c += 1,
            'G' => g += 1,
            'T' => t += 1,
            _ => {
                println!("invalid char");
                std::process::exit(1);
            }
        }
    }

    let len = (a + c + g + t) as f64;

    let tm_basic = 4.0 * (g + c) as f64 + 2.0 * (a + t) as f64;

    let na = 0.05_f64;
    let gc_percent = 100.0 * (g + c) as f64 / len;
    let tm_salt = 81.5 + 16.6 * na.log10() + 0.41 * gc_percent - 600.0 / len;

    println!("tm_simple: {:.2} °C", tm_basic);
    println!("tm_salt: {:.2} °C", tm_salt);
}
