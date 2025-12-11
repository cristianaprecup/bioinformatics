import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import threading

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

def smith_waterman_score_native(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)
    
    prev_row = [0] * (m + 1)
    curr_row = [0] * (m + 1)
    max_score = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_val = match if seq1[i-1] == seq2[j-1] else mismatch
            
            diag = prev_row[j-1] + match_val
            up = prev_row[j] + gap
            left = curr_row[j-1] + gap
            
            val = max(0, diag, up, left)
            curr_row[j] = val
            
            if val > max_score:
                max_score = val
        
        prev_row = list(curr_row)
    
    return max_score

class GenomesComparisonApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Genome Similarity Visualizer (With Scoring Metrics)")
        self.root.geometry("1200x900")

        # Layout
        self.top_panel = ttk.Frame(root, padding="10")
        self.top_panel.pack(fill="x")
        
        self.center_panel = ttk.Frame(root, padding="10")
        self.center_panel.pack(fill="both", expand=True)

        self.bottom_panel = ttk.Frame(root, padding="15", relief="sunken")
        self.bottom_panel.pack(fill="x")

        self._init_ui()

    def _init_ui(self):
        input_frame = ttk.LabelFrame(self.top_panel, text="Inputs", padding=10)
        input_frame.pack(side="left", fill="both", expand=True, padx=5)

        ttk.Label(input_frame, text="File 1 (Flu):").grid(row=0, column=0, sticky="w")
        self.file1_entry = ttk.Entry(input_frame, width=30)
        self.file1_entry.grid(row=0, column=1, padx=5)
        self.file1_entry.insert(0, "flu1.fasta")

        ttk.Label(input_frame, text="File 2 (Covid):").grid(row=1, column=0, sticky="w")
        self.file2_entry = ttk.Entry(input_frame, width=30)
        self.file2_entry.grid(row=1, column=1, padx=5)
        self.file2_entry.insert(0, "covid1.fasta")

        param_frame = ttk.LabelFrame(self.top_panel, text="Settings", padding=10)
        param_frame.pack(side="left", fill="both", expand=True, padx=5)

        ttk.Label(param_frame, text="Chunk Size:").grid(row=0, column=0)
        self.chunk_entry = ttk.Entry(param_frame, width=8)
        self.chunk_entry.insert(0, "500") 
        self.chunk_entry.grid(row=0, column=1, padx=5)
        
        control_frame = ttk.Frame(self.top_panel, padding=10)
        control_frame.pack(side="left", fill="y")
        
        self.run_btn = ttk.Button(control_frame, text="RUN COMPARISON", command=self.start_thread)
        self.run_btn.pack(fill="x", pady=5)
        
        self.progress = ttk.Progressbar(control_frame, orient="horizontal", length=200, mode="determinate")
        self.progress.pack(fill="x")

        self.status_lbl = ttk.Label(control_frame, text="Ready")
        self.status_lbl.pack()

        self.fig, self.ax = plt.subplots(figsize=(10, 5.5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.center_panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)

        ttk.Label(self.bottom_panel, text="Calculated Similarity Metrics", font=("Arial", 11, "bold")).pack(anchor="w")
        
        self.stats_text = tk.StringVar()
        self.stats_text.set("Run the analysis to see scores.")
        
        stats_lbl = ttk.Label(self.bottom_panel, textvariable=self.stats_text, font=("Consolas", 12))
        stats_lbl.pack(pady=5)

    def start_thread(self):
        self.run_btn.config(state="disabled")
        self.status_lbl.config(text="Processing...")
        
        thread = threading.Thread(target=self.run_analysis)
        thread.daemon = True
        thread.start()

    def run_analysis(self):
        f1_path = self.file1_entry.get()
        f2_path = self.file2_entry.get()
        
        seq1 = read_fasta(f1_path)
        seq2 = read_fasta(f2_path)

        if not seq1 or not seq2:
            self.root.after(0, lambda: messagebox.showerror("Error", "Files not found."))
            self.root.after(0, lambda: self.run_btn.config(state="normal"))
            return

        try:
            chunk_size = int(self.chunk_entry.get())
        except ValueError:
            return

        chunks1 = [seq1[i:i+chunk_size] for i in range(0, len(seq1), chunk_size)]
        chunks2 = [seq2[i:i+chunk_size] for i in range(0, len(seq2), chunk_size)]
        
        n = len(chunks1)
        m = len(chunks2)
        sim_matrix = np.zeros((n, m))

        for i in range(n):
            for j in range(m):
                score = smith_waterman_score_native(chunks1[i], chunks2[j])
                sim_matrix[i][j] = score
            
            prog = (i / n) * 100
            self.root.after(0, lambda v=prog: self.progress.configure(value=v))
            self.root.after(0, lambda v=int(prog): self.status_lbl.configure(text=f"{v}%"))
        
        # Eq 1: Peak Local Similarity (Max Score)
        s_max = np.max(sim_matrix)
        
        # Eq 2: Global Average Similarity
        s_avg = np.mean(sim_matrix)
        
        # Eq 3: Best Region Identity % 
        theoretical_max = chunk_size * 1 
        s_perc = (s_max / theoretical_max) * 100
        
        self.root.after(0, lambda: self.update_results(sim_matrix, chunk_size, len(seq1), len(seq2), s_max, s_avg, s_perc))

    def update_results(self, matrix, chunk_size, l1, l2, s_max, s_avg, s_perc):
        self.ax.clear()
        self.ax.set_title(f"Genome Local Similarity Heatmap\nFlu ({l1}bp) vs Covid ({l2}bp)")
        self.ax.imshow(matrix, cmap='inferno', aspect='auto', interpolation='nearest')
        self.ax.set_xlabel(f"Covid Chunks ({chunk_size}bp)")
        self.ax.set_ylabel(f"Flu Chunks ({chunk_size}bp)")
        self.canvas.draw()
        
        result_string = (
            f"1. Peak Similarity (S_max):    {int(s_max)}\n"
            f"2. Global Average (S_avg):     {s_avg:.2f}\n"
            f"3. Best Region Identity:       {s_perc:.1f}%"
        )
        self.stats_text.set(result_string)
        
        self.status_lbl.config(text="Done!")
        self.run_btn.config(state="normal")
        self.progress['value'] = 100

if __name__ == "__main__":
    root = tk.Tk()
    app = GenomesComparisonApp(root)
    root.mainloop()