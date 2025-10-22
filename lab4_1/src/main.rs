use std::io::{self, Write};

fn main() {
    print!("Enter RNA sequence: ");
    io::stdout().flush().unwrap(); 

    let mut input = String::new();
    io::stdin().read_line(&mut input).unwrap();

    let rna: String = input
        .trim()
        .to_uppercase()
        .replace('T', "U")
        .chars()
        .filter(|c| matches!(c, 'A' | 'C' | 'G' | 'U'))
        .collect();

    let Some(mut i) = rna.find("AUG") else {
        println!("No start codon found.");
        return;
    };

    let bytes = rna.as_bytes();
    let mut nume: Vec<&'static str> = Vec::new();

    while i + 3 <= bytes.len() {
        let codon = std::str::from_utf8(&bytes[i..i + 3]).unwrap();

        if matches!(codon, "UAA" | "UAG" | "UGA") {
            break;
        }

        if let Some((n, l)) = codon_to_aa(codon) {
            nume.push(n);
        } else {
            println!("Invalid codon: {}", codon);
            break;
        }
        i += 3;
    }

    if nume.is_empty() {
        println!("No amino acids translated between AUG and Stop.");
    } else {
        println!("Amino acids: {}", nume.join("-"));
    }
}

fn codon_to_aa(c: &str) -> Option<(&'static str, char)> {
    Some(match c {
        "UUU" | "UUC" => ("Phe", 'F'),
        "UUA" | "UUG" | "CUU" | "CUC" | "CUA" | "CUG" => ("Leu", 'L'),
        "AUU" | "AUC" | "AUA" => ("Ile", 'I'),
        "AUG" => ("Met", 'M'),
        "GUU" | "GUC" | "GUA" | "GUG" => ("Val", 'V'),
        "UCU" | "UCC" | "UCA" | "UCG" | "AGU" | "AGC" => ("Ser", 'S'),
        "CCU" | "CCC" | "CCA" | "CCG" => ("Pro", 'P'),
        "ACU" | "ACC" | "ACA" | "ACG" => ("Thr", 'T'),
        "GCU" | "GCC" | "GCA" | "GCG" => ("Ala", 'A'),
        "UAU" | "UAC" => ("Tyr", 'Y'),
        "CAU" | "CAC" => ("His", 'H'),
        "CAA" | "CAG" => ("Gln", 'Q'),
        "AAU" | "AAC" => ("Asn", 'N'),
        "AAA" | "AAG" => ("Lys", 'K'),
        "GAU" | "GAC" => ("Asp", 'D'),
        "GAA" | "GAG" => ("Glu", 'E'),
        "UGU" | "UGC" => ("Cys", 'C'),
        "UGG" => ("Trp", 'W'),
        "CGU" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => ("Arg", 'R'),
        "GGU" | "GGC" | "GGA" | "GGG" => ("Gly", 'G'),
        _ => return None,
    })
}
