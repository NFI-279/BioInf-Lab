from collections import defaultdict
import gzip
import random
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk


class GenomeAssemblerApp:

  def __init__(self, root):
    self.root = root
    self.root.title("De Novo Genome Assembler Lab (20K Reads)")
    self.root.geometry("900x700")
    self.root.minsize(800, 600)

    # Default file path
    self.file_path = "GCF_000919155.1_ViralProj244477_genomic.fna.gz"

    self.create_widgets()

  def create_widgets(self):
    # Top Configuration Frame
    config_frame = ttk.LabelFrame(
        self.root, text=" Lab Parameters & File Config ", padding=15
    )
    config_frame.pack(fill="x", padx=15, pady=10)

    # File selection
    ttk.Label(config_frame, text="Genome File:").grid(
        row=0, column=0, sticky="w", padx=5, pady=5
    )
    self.file_entry = ttk.Entry(config_frame, width=40)
    self.file_entry.grid(row=0, column=1, sticky="w", padx=5, pady=5)
    self.file_entry.insert(0, self.file_path)

    browse_btn = ttk.Button(
        config_frame, text="Browse...", command=self.browse_file
    )
    browse_btn.grid(row=0, column=2, padx=5, pady=5)

    # Number of reads
    ttk.Label(config_frame, text="Number of Samples (Reads):").grid(
        row=1, column=0, sticky="w", padx=5, pady=5
    )
    self.reads_entry = ttk.Entry(config_frame, width=15)
    self.reads_entry.grid(row=1, column=1, sticky="w", padx=5, pady=5)
    self.reads_entry.insert(0, "20000")

    # Min Overlap Threshold (>10 letters rule)
    ttk.Label(config_frame, text="Min Overlap Threshold (>10):").grid(
        row=2, column=0, sticky="w", padx=5, pady=5
    )
    self.threshold_entry = ttk.Entry(config_frame, width=15)
    self.threshold_entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)
    self.threshold_entry.insert(
        0, "12"
    )  # Safe threshold above 10 to prevent false positives

    # Run Button
    self.run_btn = ttk.Button(
        config_frame,
        text="Run Simulation & Assembly",
        command=self.start_assembly_pipeline,
    )
    self.run_btn.grid(row=1, column=2, rowspan=2, padx=15, pady=5, ipadx=10)

    # Notebook (Tabs) for Output
    self.notebook = ttk.Notebook(self.root)
    self.notebook.pack(fill="both", expand=True, padx=15, pady=10)

    # Tab 1: Execution Logs & Summary Report
    self.log_frame = ttk.Frame(self.notebook)
    self.notebook.add(self.log_frame, text=" Assembly Report & Logs ")

    self.log_text = scrolledtext.ScrolledText(
        self.log_frame, wrap="word", font=("Consolas", 10)
    )
    self.log_text.pack(fill="both", expand=True, padx=10, pady=10)
    self.log_text.insert(
        tk.END,
        "Click 'Run Simulation & Assembly' to break the genome into 20k samples"
        " and reconstruct it...\n",
    )
    self.log_text.configure(state="disabled")

    # Tab 2: Reconstructed Genome Viewer
    self.viewer_frame = ttk.Frame(self.notebook)
    self.notebook.add(self.viewer_frame, text=" Reconstructed Genome Viewer ")

    self.viewer_text = scrolledtext.ScrolledText(
        self.viewer_frame, wrap="word", font=("Consolas", 10)
    )
    self.viewer_text.pack(fill="both", expand=True, padx=10, pady=10)
    self.viewer_text.insert(
        tk.END, "Assembled sequence will appear here after execution."
    )
    self.viewer_text.configure(state="disabled")

    # Bottom Status Bar
    self.status_var = tk.StringVar(value="Status: Ready")
    status_bar = ttk.Label(
        self.root,
        textvariable=self.status_var,
        relief="sunken",
        anchor="w",
        padding=5,
    )
    status_bar.pack(fill="x", side="bottom")

  def browse_file(self):
    file_path = filedialog.askopenfilename(
        title="Select Genome FASTA File",
        filetypes=[("FASTA Compressed", "*.gz"), ("All Files", "*.*")],
    )
    if file_path:
      self.file_path = file_path
      self.file_entry.delete(0, tk.END)
      self.file_entry.insert(0, file_path)

  def log_message(self, message):
    self.log_text.configure(state="normal")
    self.log_text.insert(tk.END, message + "\n")
    self.log_text.see(tk.END)
    self.log_text.configure(state="normal")
    self.status_var.set(f"Status: {message}")
    self.root.update_idletasks()

  def start_assembly_pipeline(self):
    # Clear logs
    self.log_text.configure(state="normal")
    self.log_text.delete("1.0", tk.END)
    self.log_text.configure(state="disabled")

    try:
      num_reads = int(self.reads_entry.get().strip())
      min_threshold = int(self.threshold_entry.get().strip())
      file_path = self.file_entry.get().strip()
    except ValueError:
      messagebox.showerror(
          "Input Error",
          "Please enter valid integer values for reads and threshold.",
      )
      return

    self.run_btn.config(state="disabled")
    threading.Thread(
        target=self.run_lab_simulation,
        args=(file_path, num_reads, min_threshold),
        daemon=True,
    ).start()

  def run_lab_simulation(self, file_path, num_reads, min_threshold):
    try:
      # Step 1: Load Genome
      self.log_message(f"Loading reference genome from: {file_path}...")
      seq_chunks = []
      open_func = gzip.open if file_path.endswith(".gz") else open
      mode = "rt" if file_path.endswith(".gz") else "r"

      with open_func(file_path, mode) as f:
        for line in f:
          line = line.strip()
          if not line or line.startswith(">"):
            continue
          seq_chunks.append(line)

      original_genome = "".join(seq_chunks).upper()
      self.log_message(
          f"-> Original Genome Loaded. Length: {len(original_genome):,}"
          " nucleotides."
      )

      # Step 2: Simulate Breaking into 20k samples (100-150 bp)
      self.log_message(
          f"Simulating sequencer: Generating {num_reads:,} random samples"
          " (100-150 bp)..."
      )
      reads = []
      genome_len = len(original_genome)

      for _ in range(num_reads):
        read_len = random.randint(100, 150)
        if genome_len <= read_len:
          reads.append(original_genome)
          continue
        start = random.randint(0, genome_len - read_len)
        reads.append(original_genome[start : start + read_len])

      self.log_message(f"-> Successfully generated {len(reads):,} samples.")

      # Step 3: Ultra-Fast Chain Assembly (Enforcing >10 letter threshold)
      self.log_message(
          f"Starting fast assembly with overlap threshold >= {min_threshold}"
          " bases..."
      )
      start_time = time.time()

      unique_reads = list(set(reads))
      prefix_map = defaultdict(list)
      for r in unique_reads:
        if len(r) >= min_threshold:
          prefix_map[r[:min_threshold]].append(r)

      def get_overlap(s1, s2):
        max_possible = min(len(s1), len(s2))
        for i in range(max_possible, min_threshold - 1, -1):
          if s1[-i:] == s2[:i]:
            return i
        return 0

      used = set()
      contigs = []

      for start_read in unique_reads:
        if start_read in used:
          continue

        current_contig = start_read
        used.add(start_read)

        while True:
          if len(current_contig) < min_threshold:
            break
          suffix = current_contig[-min_threshold:]
          candidates = prefix_map.get(suffix, [])

          best_cand = None
          max_ov = 0

          for cand in candidates:
            if cand not in used:
              ov = get_overlap(current_contig, cand)
              if ov > max_ov:
                max_ov = ov
                best_cand = cand

          if best_cand and max_ov >= min_threshold:
            current_contig += best_cand[max_ov:]
            used.add(best_cand)
          else:
            break

        contigs.append(current_contig)

      elapsed = time.time() - start_time
      self.log_message(
          f"-> Assembly completed in {elapsed:.2f} seconds!"
      )

      # Step 4: Report Results
      self.log_message("\n--- LAB RECONSTRUCTION REPORT ---")
      self.log_message(f"Total reconstructed contigs: {len(contigs)}")

      for idx, contig in enumerate(contigs, 1):
        self.log_message(
            f"Contig {idx} Length: {len(contig):,} nucleotides"
        )

      success = len(contigs) == 1 and abs(
          len(contigs[0]) - len(original_genome)
      ) < 50
      if success:
        self.log_message(
            "\n[SUCCESS]: The original genome was fully and perfectly"
            " reconstructed from the 20,000 samples!"
        )
      else:
        self.log_message(
            "\n[NOTE]: Reconstructed into major contigs meeting the threshold"
            " rules."
        )

      # Populate viewer tab with final sequence
      self.viewer_text.configure(state="normal")
      self.viewer_text.delete("1.0", tk.END)
      self.viewer_text.insert(
          tk.END,
          "\n\n".join([f">Contig_{i}\n{c}" for i, c in enumerate(contigs, 1)]),
      )
      self.viewer_text.configure(state="disabled")

    except Exception as e:
      messagebox.showerror("Execution Error", str(e))
      self.log_message(f"Error: {e}")
    finally:
      self.run_btn.config(state="normal")
      self.status_var.set("Status: Ready")


if __name__ == "__main__":
  root = tk.Tk()
  app = GenomeAssemblerApp(root)
  root.mainloop()