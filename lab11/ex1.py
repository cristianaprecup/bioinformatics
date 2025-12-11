import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import ListedColormap

class DNAAlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Needleman-Wunsch Alignment Tool")
        self.root.geometry("1100x700")

        self.left_panel = ttk.Frame(root, padding="10")
        self.left_panel.grid(row=0, column=0, sticky="nsew")

        self.right_panel = ttk.Frame(root, padding="10")
        self.right_panel.grid(row=0, column=1, sticky="nsew")
        
        self.bottom_panel = ttk.Frame(root, padding="10")
        self.bottom_panel.grid(row=1, column=0, columnspan=2, sticky="nsew")

        self._init_inputs()
        self._init_graphs()
        self._init_output()
        
        self.run_alignment()

    def _init_inputs(self):
        seq_frame = ttk.LabelFrame(self.left_panel, text="Sequences", padding=10)
        seq_frame.pack(fill="x", pady=5)

        ttk.Label(seq_frame, text="Sq 1:").pack(anchor="w")
        self.seq1_entry = ttk.Entry(seq_frame, width=30)
        self.seq1_entry.pack(fill="x", pady=(0, 5))
        self.seq1_entry.insert(0, "ACCGTGAAGCCAATAC")

        ttk.Label(seq_frame, text="Sq 2:").pack(anchor="w")
        self.seq2_entry = ttk.Entry(seq_frame, width=30)
        self.seq2_entry.pack(fill="x")
        self.seq2_entry.insert(0, "AGCGTGCAGCCAATAC")

        param_frame = ttk.LabelFrame(self.left_panel, text="Parameters", padding=10)
        param_frame.pack(fill="x", pady=10)

        ttk.Label(param_frame, text="Gap:").grid(row=0, column=0, sticky="e")
        self.gap_entry = ttk.Entry(param_frame, width=5)
        self.gap_entry.grid(row=0, column=1, padx=5, pady=2)
        self.gap_entry.insert(0, "0")

        ttk.Label(param_frame, text="Match:").grid(row=1, column=0, sticky="e")
        self.match_entry = ttk.Entry(param_frame, width=5)
        self.match_entry.grid(row=1, column=1, padx=5, pady=2)
        self.match_entry.insert(0, "1")

        ttk.Label(param_frame, text="Missmatch:").grid(row=2, column=0, sticky="e")
        self.mismatch_entry = ttk.Entry(param_frame, width=5)
        self.mismatch_entry.grid(row=2, column=1, padx=5, pady=2)
        self.mismatch_entry.insert(0, "-1")

        self.align_btn = ttk.Button(self.left_panel, text="Align", command=self.run_alignment)
        self.align_btn.pack(fill="x", pady=20)

    def _init_graphs(self):
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.right_panel)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.fig.tight_layout()

    def _init_output(self):
        ttk.Label(self.bottom_panel, text="Alignment Output:").pack(anchor="w")
        self.output_text = scrolledtext.ScrolledText(self.bottom_panel, height=10, font=("Courier", 10))
        self.output_text.pack(fill="both", expand=True)

    def run_alignment(self):
        s1 = self.seq1_entry.get() 
        s2 = self.seq2_entry.get() 
        
        try:
            gap = int(self.gap_entry.get())
            match = int(self.match_entry.get())
            mismatch = int(self.mismatch_entry.get())
        except ValueError:
            return

        n_cols = len(s1)
        m_rows = len(s2)
        
        score_matrix = np.zeros((m_rows + 1, n_cols + 1), dtype=int)
        
        for i in range(m_rows + 1): score_matrix[i][0] = i * gap
        for j in range(n_cols + 1): score_matrix[0][j] = j * gap

        for i in range(1, m_rows + 1):
            for j in range(1, n_cols + 1):
                diag_score = score_matrix[i-1][j-1] + (match if s2[i-1] == s1[j-1] else mismatch)
                up_score = score_matrix[i-1][j] + gap    
                left_score = score_matrix[i][j-1] + gap  
                
                score_matrix[i][j] = max(diag_score, up_score, left_score)

        align1, align2 = "", ""
        
        traceback_visual = np.zeros((m_rows + 1, n_cols + 1))
        
        i, j = m_rows, n_cols
        traceback_visual[i][j] = 1 

        while i > 0 or j > 0:
            current = score_matrix[i][j]
            
            score_diag = -9999
            if i > 0 and j > 0:
                val = match if s2[i-1] == s1[j-1] else mismatch
                score_diag = score_matrix[i-1][j-1] + val
            
            score_up = score_matrix[i-1][j] + gap if i > 0 else -9999
            
            if i > 0 and j > 0 and current == score_diag:
                align1 = s1[j-1] + align1
                align2 = s2[i-1] + align2
                i -= 1; j -= 1
            elif i > 0 and current == score_up:
                align1 = "-" + align1
                align2 = s2[i-1] + align2
                i -= 1
            else:
                align1 = s1[j-1] + align1
                align2 = "-" + align2
                j -= 1
            
            traceback_visual[i][j] = 1 

        matches = sum(1 for a, b in zip(align1, align2) if a == b)
        similarity = (matches / len(align1)) * 100
        bars = "".join(["|" if a == b else " " for a, b in zip(align1, align2)])

        txt = f"Show Alignment:\n\n{align1}\n{bars}\n{align2}\n\n"
        txt += f"Matches = {matches}\nLength = {len(align1)}\nSimilarity = {int(similarity)} %"
        
        self.output_text.delete(1.0, tk.END)
        self.output_text.insert(tk.END, txt)

        self.ax1.clear()
        self.ax1.set_title("Score Matrix")
        self.ax1.imshow(score_matrix, cmap='magma', aspect='auto')
        self.ax1.set_xlabel(f"Seq 1 ({len(s1)})")
        self.ax1.set_ylabel(f"Seq 2 ({len(s2)})")

        self.ax2.clear()
        self.ax2.set_title("Traceback Path")
        
        cmap_trace = ListedColormap(['#FFFFE0', '#D32F2F']) 
        
        self.ax2.imshow(traceback_visual, cmap=cmap_trace, aspect='auto', interpolation='nearest')
        
        self.ax2.set_xticks(np.arange(-.5, n_cols + 1, 1), minor=True)
        self.ax2.set_yticks(np.arange(-.5, m_rows + 1, 1), minor=True)
        self.ax2.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
        self.ax2.tick_params(which='minor', size=0)
        
        self.ax2.set_xlabel("Seq 1")
        
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = DNAAlignmentApp(root)
    root.mainloop()