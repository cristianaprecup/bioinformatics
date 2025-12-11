import tkinter as tk
from tkinter import ttk, messagebox
from tkinter import scrolledtext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os

def read_fasta(file_path):
    if not os.path.exists(file_path):
        return None
    
    sequence = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if not line.startswith(">"):
                sequence.append(line)
    return "".join(sequence).upper()

def smith_waterman_score_only(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)
    
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    max_score = 0
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_val = match if seq1[i-1] == seq2[j-1] else mismatch
            
            diag = score_matrix[i-1][j-1] + match_val
            up = score_matrix[i-1][j] + gap
            left = score_matrix[i][j-1] + gap
            
            cell_score = max(0, diag, up, left)
            score_matrix[i][j] = cell_score
            
            if cell_score > max_score:
                max_score = cell_score
                
    return max_score

class GenomesComparisonApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Similarity Visualizer (Local Alignment)")
        self.root.geometry("1200x800")

        self.top_panel = ttk.Frame(root, padding="10")
        self.top_panel.pack(fill="x")
        
        self.center_panel = ttk.Frame(root, padding="10")
        self.center_panel.pack(fill="both", expand=True)

        self._init_ui()

    def _init_ui(self):
        input_frame = ttk.LabelFrame(self.top_panel, text="Genome Inputs", padding=10)
        input_frame.pack(side="left", fill="both", expand=True, padx=5)

        ttk.Label(input_frame, text="File 1 (Flu):").grid(row=0, column=0, sticky="w")
        self.file1_entry = ttk.Entry(input_frame, width=30)
        self.file1_entry.grid(row=0, column=1, padx=5)
        self.file1_entry.insert(0, "flu1.fasta")

        ttk.Label(input_frame, text="File 2 (Covid):").grid(row=1, column=0, sticky="w")
        self.file2_entry = ttk.Entry(input_frame, width=30)
        self.file2_entry.grid(row=1, column=1, padx=5)
        self.file2_entry.insert(0, "covid1.fasta")

        param_frame = ttk.LabelFrame(self.top_panel, text="Analysis Parameters", padding=10)
        param_frame.pack(side="left", fill="both", expand=True, padx=5)

        ttk.Label(param_frame, text="Chunk Size (bp):").grid(row=0, column=0, sticky="w")
        self.chunk_entry = ttk.Entry(param_frame, width=8)
        self.chunk_entry.insert(0, "100") 
        self.chunk_entry.grid(row=0, column=1)
        
        ttk.Label(param_frame, text="(Smaller = Higher Res, Slower)").grid(row=0, column=2, sticky="w", padx=5)

        btn_frame = ttk.Frame(self.top_panel, padding=10)
        btn_frame.pack(side="left", fill="y")
        self.run_btn = ttk.Button(btn_frame, text="RUN COMPARISON", command=self.run_analysis)
        self.run_btn.pack(fill="both", expand=True)

        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.center_panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

    def run_analysis(self):
        f1_path = self.file1_entry.get()
        f2_path = self.file2_entry.get()
        
        seq1 = read_fasta(f1_path)
        seq2 = read_fasta(f2_path)

        if not seq1 or not seq2:
            messagebox.showerror("Error", "Could not read FASTA files.\nMake sure ./flu1.fasta and ./covid1.fasta exist.")
            return

        try:
            chunk_size = int(self.chunk_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Chunk size must be an integer.")
            return

        chunks1 = [seq1[i:i+chunk_size] for i in range(0, len(seq1), chunk_size)]
        chunks2 = [seq2[i:i+chunk_size] for i in range(0, len(seq2), chunk_size)]

        n_chunks = len(chunks1)
        m_chunks = len(chunks2)
        
        if n_chunks * m_chunks > 100000:
             if not messagebox.askyesno("Warning", f"This will generate {n_chunks*m_chunks} comparisons. It might be slow. Continue?"):
                 return

        sim_matrix = np.zeros((n_chunks, m_chunks))

        print(f"Starting comparison: {n_chunks} x {m_chunks} blocks.")
        
        for i in range(n_chunks):
            for j in range(m_chunks):
                s1_segment = chunks1[i]
                s2_segment = chunks2[j]
                
                raw_score = smith_waterman_score_only(s1_segment, s2_segment)
                sim_matrix[i][j] = raw_score
            
            if i % 10 == 0:
                self.root.update() 

        self.ax.clear()
        self.ax.set_title(f"Genome Similarity Map (Local Alignment Scores)\nFlu ({len(seq1)}bp) vs Covid ({len(seq2)}bp)")
        
        cax = self.ax.imshow(sim_matrix, cmap='inferno', aspect='auto', interpolation='nearest')
        
        self.ax.set_xlabel(f"Covid-19 Chunks (Size {chunk_size})")
        self.ax.set_ylabel(f"Influenza Chunks (Size {chunk_size})")
        
        if not hasattr(self, 'cbar'):
            self.cbar = self.fig.colorbar(cax, ax=self.ax)
        else:
            self.cbar.update_normal(cax)

        self.canvas.draw()
        messagebox.showinfo("Done", "Comparison Complete.\nBrighter regions indicate higher genetic similarity.")

if __name__ == "__main__":
    root = tk.Tk()
    app = GenomesComparisonApp(root)
    root.mainloop() 