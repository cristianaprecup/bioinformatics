#make an application that is able to find the alphabet of a sequence of text. this seq may be an ARN or ADN or protein seq

def find_alphabet(sequence):
    sequence = sequence.upper()
    rna_alphabet = {'A', 'U', 'G', 'C'}
    dna_alphabet = {'A', 'T', 'G', 'C'}
    protein_alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    if set(sequence).issubset(rna_alphabet):
        return "RNA sequence"
    elif set(sequence).issubset(dna_alphabet):
        return "DNA sequence"
    elif set(sequence).issubset(protein_alphabet):
        return "Protein sequence"
    else:
        return "Unknown sequence type"


if __name__ == "__main__":
    seq = input("Enter a sequence: ").strip()
    result = find_alphabet(seq)
    print(f"The sequence is identified as: {result}")