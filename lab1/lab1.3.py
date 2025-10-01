import tkinter as tk
from tkinter import filedialog, messagebox

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def calculate_percentage(sequence):
    total_length = len(sequence)
    letter_count = {}

    for letter in sequence:
        if letter in letter_count:
            letter_count[letter] += 1
        else:
            letter_count[letter] = 1

    letter_percentage = {
        letter: (count / total_length) * 100
        for letter, count in letter_count.items()
    }

    return letter_percentage

def open_fasta_file():
    file_path = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=(("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
    )
    if not file_path:
        return

    try:
        sequence = read_fasta(file_path)
        percentages = calculate_percentage(sequence)

        result_text.delete(1.0, tk.END)
        result_text.insert(tk.END, "Letter percentages:\n")
        for letter, percentage in percentages.items():
            result_text.insert(tk.END, f"{letter}: {percentage:.2f}%\n")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {e}")

def main():
    global result_text
    root = tk.Tk()
    root.title("FASTA File Analyzer")

    open_button = tk.Button(root, text="Open FASTA File", command=open_fasta_file)
    open_button.pack(pady=10)

    result_text = tk.Text(root, width=50, height=20)
    result_text.pack(pady=10)

    root.mainloop()

if __name__ == "__main__":
    main()