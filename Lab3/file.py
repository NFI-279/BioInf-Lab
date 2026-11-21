from collections import Counter
import gzip
import itertools
import os
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from Bio import Entrez, SeqIO
import openai

# Configure NCBI Entrez (Required by NCBI)
Entrez.email = "your_email@example.com"


def generate_rna_codons():
  """Generates all 64 RNA codon combinations and initializes dictionary."""
  nucleotides = ["A", "C", "G", "U"]
  all_codons = ["".join(p) for p in itertools.product(nucleotides, repeat=3)]
  return {codon: 0 for codon in all_codons}


class VirusDietApp:

  def __init__(self, root):
    self.root = root
    self.root.title("Viral Codon Genome & Dietary Advisor")
    self.root.geometry("900x750")
    self.root.minsize(800, 650)

    self.local_file_path = None
    self.create_widgets()

  def create_widgets(self):
    # Top Frame: Inputs
    input_frame = ttk.LabelFrame(
        self.root, text=" 1. Target Selection & API Config ", padding=15
    )
    input_frame.pack(fill="x", padx=15, pady=10)

    ttk.Label(input_frame, text="Virus Name / Taxon:").grid(
        row=0, column=0, sticky="w", padx=5, pady=5
    )
    self.virus_entry = ttk.Entry(input_frame, width=25)
    self.virus_entry.grid(row=0, column=1, sticky="w", padx=5, pady=5)
    self.virus_entry.insert(0, "Influenza A virus")

    # Local File Selection Option
    self.file_label = ttk.Label(
        input_frame, text="Or Use Local File: None selected", foreground="gray"
    )
    self.file_label.grid(row=1, column=0, columnspan=2, sticky="w", padx=5, pady=2)

    browse_btn = ttk.Button(
        input_frame, text="Browse .fasta.gz", command=self.browse_file
    )
    browse_btn.grid(row=1, column=2, sticky="w", padx=5, pady=2)

    ttk.Label(input_frame, text="OpenAI API Key (Optional):").grid(
        row=2, column=0, sticky="w", padx=5, pady=5
    )
    self.api_key_entry = ttk.Entry(input_frame, width=25, show="*")
    self.api_key_entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)

    self.run_button = ttk.Button(
        input_frame, text="Run Analysis", command=self.start_analysis
    )
    self.run_button.grid(row=0, column=2, rowspan=2, padx=15, pady=5, ipadx=10)

    # Notebook (Tabs) for Results
    self.notebook = ttk.Notebook(self.root)
    self.notebook.pack(fill="both", expand=True, padx=15, pady=10)

    # Tab 1: Top 10 Codons Table
    self.codon_frame = ttk.Frame(self.notebook)
    self.notebook.add(self.codon_frame, text=" Top 10 Codon Frequencies ")

    self.tree = ttk.Treeview(
        self.codon_frame, columns=("Rank", "Codon", "Frequency"), show="headings"
    )
    self.tree.heading("Rank", text="Rank")
    self.tree.heading("Codon", text="Codon")
    self.tree.heading("Frequency", text="Frequency Count")
    self.tree.column("Rank", width=80, anchor="center")
    self.tree.column("Codon", width=150, anchor="center")
    self.tree.column("Frequency", width=200, anchor="center")

    tree_scroll = ttk.Scrollbar(
        self.codon_frame, orient="vertical", command=self.tree.yview
    )
    self.tree.configure(yscrollcommand=tree_scroll.set)
    self.tree.pack(side="left", fill="both", expand=True, padx=5, pady=5)
    tree_scroll.pack(side="right", fill="y", pady=5)

    # Tab 2: LLM Dietary Recommendations
    self.diet_frame = ttk.Frame(self.notebook)
    self.notebook.add(
        self.diet_frame, text=" AI Dietary Advice (Eat vs. Avoid) "
    )

    self.diet_text = scrolledtext.ScrolledText(
        self.diet_frame, wrap="word", font=("Arial", 11)
    )
    self.diet_text.pack(fill="both", expand=True, padx=10, pady=10)
    self.diet_text.insert(
        tk.END,
        "Run an analysis to generate tailored dietary recommendations based"
        " on genomic composition.",
    )
    self.diet_text.configure(state="disabled")

    # Bottom Frame: Status Bar
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
        title="Select Genomic FASTA File",
        filetypes=[
            ("FASTA / Compressed Files", "*.fna *.fasta *.fa *.gz"),
            ("All Files", "*.*"),
        ],
    )
    if file_path:
      self.local_file_path = file_path
      filename_only = os.path.basename(file_path)
      self.file_label.config(
          text=f"Using Local: {filename_only}", foreground="green"
      )

  def log_status(self, message):
    self.status_var.set(f"Status: {message}")
    self.root.update_idletasks()

  def start_analysis(self):
    virus_name = self.virus_entry.get().strip()
    if not self.local_file_path and not virus_name:
      messagebox.showerror(
          "Error", "Please enter a virus name or select a local file."
      )
      return

    self.run_button.config(state="disabled")
    threading.Thread(
        target=self.run_pipeline, args=(virus_name,), daemon=True
    ).start()

  def run_pipeline(self, virus_name):
    try:
      fasta_filename = None

      # Determine if we use local file or fetch from NCBI
      if self.local_file_path:
        self.log_status(
            f"Reading local file: {os.path.basename(self.local_file_path)}..."
        )
        fasta_filename = self.local_file_path
        target_name = (
            virus_name if virus_name else os.path.basename(self.local_file_path)
        )
      else:
        self.log_status(
            f"Searching NCBI Entrez for '{virus_name}' reference genome..."
        )
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{virus_name} complete genome[Title]",
            retmax=1,
        )
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
          handle = Entrez.esearch(db="nucleotide", term=virus_name, retmax=1)
          record = Entrez.read(handle)
          handle.close()

        if not record["IdList"]:
          raise ValueError(
              f"No genomic records found on NCBI for '{virus_name}'."
          )

        ncbi_id = record["IdList"][0]
        self.log_status(f"Fetching FASTA sequence for NCBI ID: {ncbi_id}...")

        fetch_handle = Entrez.efetch(
            db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text"
        )
        fasta_data = fetch_handle.read()
        fetch_handle.close()

        fasta_filename = f"{ncbi_id}_genomic.fna"
        with open(fasta_filename, "w") as f:
          f.write(fasta_data)
        target_name = virus_name

      # Parse sequence and analyze codons (handling both normal and .gz files)
      self.log_status("Running sliding window codon frequency analysis...")
      codon_counts = generate_rna_codons()
      seq_chunks = []

      open_func = gzip.open if fasta_filename.endswith(".gz") else open
      mode = "rt" if fasta_filename.endswith(".gz") else "r"

      with open_func(fasta_filename, mode) as f:
        for line in f:
          line = line.strip()
          if not line or line.startswith(">"):
            continue
          seq_chunks.append(line)

      full_seq = "".join(seq_chunks).upper().replace("T", "U")

      # Sliding window of 3
      valid_nucs = {"A", "C", "G", "U"}
      for i in range(len(full_seq) - 2):
        codon = full_seq[i : i + 3]
        if all(n in valid_nucs for n in codon):
          if codon in codon_counts:
            codon_counts[codon] += 1

      top_10 = sorted(
          codon_counts.items(), key=lambda item: item[1], reverse=True
      )[:10]

      self.root.after(0, lambda: self.update_table(top_10))

      # Query LLM for Dietary Guidance
      self.log_status("Generating AI dietary recommendations...")
      api_key = self.api_key_entry.get().strip()
      diet_advice = self.get_llm_recommendations(target_name, top_10, api_key)

      self.root.after(0, lambda: self.update_diet_tab(diet_advice))
      self.log_status("Analysis complete successfully!")

    except Exception as e:
      self.root.after(
          0, lambda: messagebox.showerror("Execution Error", str(e))
      )
      self.log_status("Error encountered during execution.")
    finally:
      self.root.after(0, lambda: self.run_button.config(state="normal"))

  def update_table(self, top_10):
    for row in self.tree.get_children():
      self.tree.delete(row)
    for rank, (codon, freq) in enumerate(top_10, start=1):
      self.tree.insert("", "end", values=(rank, codon, f"{freq:,}"))

  def get_llm_recommendations(self, target_name, top_10, api_key):
    prompt = f"""
    You are an expert clinical nutritionist and molecular biologist. 
    A patient is researching dietary mitigation for: {target_name}.
    We extracted the genomic sequence and found the top 10 most frequent RNA codons:
    {top_10}
    
    Based on viral replication mechanics and amino acid demands:
    1. Explain what foods/nutrients contain amino acids or support biological pathways corresponding to these high-frequency codons.
    2. Provide a clear breakdown of what to **EAT** to boost immune response.
    3. Provide a clear breakdown of what to **AVOID**.
    
    Format the response clearly with bold headings and bullet points.
    """

    if api_key:
      try:
        client = openai.OpenAI(api_key=api_key)
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{
                "role": "user",
                "content": prompt,
            }],
        )
        return response.choices[0].message.content
      except Exception as e:
        return (
            f"OpenAI API call failed: {e}\n\nShowing heuristic nutritional"
            " template instead."
        )

    return f"""
    [Simulated AI Nutritional Report for: {target_name}]
    
    Top Codons Analyzed: {', '.join([c[0] for c in top_10[:5]])}
    
    ## 🥗 What to EAT (Immune Defense & Cellular Support)
    * **High-Antioxidant Produce:** Berries, leafy greens, and citrus fruits to combat oxidative stress.
    * **Lean Proteins & Essential Amino Acids:** Chicken bone broth, wild-caught fish, and eggs for cellular repair.
    * **Anti-Inflammatory Fats:** Omega-3 fatty acids from salmon, flaxseeds, and walnuts.
    * **Micronutrients (Zinc & Vitamin D):** Pumpkin seeds and mushrooms.

    ## 🚫 What to AVOID
    * **Refined Sugars & Processed Foods:** Temporarily impairs immune macrophage activity.
    * **Excessive Alcohol:** Suppresses mucosal immunity.
    * **Pro-inflammatory Seed Oils:** Can aggravate inflammatory pathways.
    """

  def update_diet_tab(self, advice):
    self.diet_text.configure(state="normal")
    self.diet_text.delete("1.0", tk.END)
    self.diet_text.insert(tk.END, advice)
    self.diet_text.configure(state="disabled")


if __name__ == "__main__":
  root = tk.Tk()
  app = VirusDietApp(root)
  root.mainloop()