"""Restriction-enzyme digestion lab GUI.

Run with:  python restriction_digest_gui.py
The program accepts ordinary NCBI FASTA files and reports all seven
non-empty combinations of EcoRI, BamHI, and HindIII.
"""
from __future__ import annotations

import itertools
import gzip
import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


ENZYMES = {
    "EcoRI": ("GAATTC", 1),       # G|AATTC
    "BamHI": ("GGATCC", 1),       # G|GATCC
    "HindIII": ("AAGCTT", 1),     # A|AGCTT
}


def read_fasta(path: str) -> tuple[str, str]:
    """Return (description, uppercase sequence), ignoring whitespace."""
    source = Path(path)
    if source.suffix.lower() == ".gz":
        with gzip.open(source, "rt", encoding="utf-8-sig") as handle:
            text = handle.read()
    else:
        text = source.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    headers = [line[1:].strip() for line in lines if line.startswith(">")]
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    sequence = re.sub(r"\s+", "", sequence).upper().replace("U", "T")
    if not sequence or re.search(r"[^ACGTN]", sequence):
        raise ValueError("The FASTA sequence must contain only A, C, G, T, or N.")
    return (headers[0] if headers else Path(path).name), sequence


def find_cuts(sequence: str, enzyme: str) -> list[int]:
    pattern, offset = ENZYMES[enzyme]
    # A cut coordinate is 1-based and means: cut between coordinate and coordinate+1.
    return [match.start() + offset for match in re.finditer(f"(?={pattern})", sequence)]


def digest(sequence: str, enzymes: tuple[str, ...]) -> tuple[list[int], list[int]]:
    cuts = sorted({cut for enzyme in enzymes for cut in find_cuts(sequence, enzyme)})
    boundaries = [0, *cuts, len(sequence)]
    fragments = [boundaries[i + 1] - boundaries[i] for i in range(len(boundaries) - 1)]
    return cuts, fragments


class DigestApp(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Restriction Enzyme Digestion Lab")
        self.geometry("1100x720")
        self.minsize(850, 580)
        self.sequence = ""
        self.description = ""
        self.path_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Choose an NCBI FASTA file to begin.")
        self._build_ui()

    def _build_ui(self) -> None:
        top = ttk.LabelFrame(self, text="1. Input FASTA genome")
        top.pack(fill="x", padx=10, pady=(10, 5))
        ttk.Entry(top, textvariable=self.path_var).pack(side="left", fill="x", expand=True, padx=8, pady=8)
        ttk.Button(top, text="Browse…", command=self.browse).pack(side="left", padx=4)
        ttk.Button(top, text="Load FASTA", command=self.load).pack(side="left", padx=(0, 8))

        controls = ttk.LabelFrame(self, text="2. Enzymes")
        controls.pack(fill="x", padx=10, pady=5)
        self.vars = {name: tk.BooleanVar(value=True) for name in ENZYMES}
        for name, (site, offset) in ENZYMES.items():
            ttk.Checkbutton(controls, text=f"{name}  {site}  (cut after base {offset})",
                            variable=self.vars[name]).pack(side="left", padx=10, pady=8)
        ttk.Button(controls, text="Analyze all 7 combinations", command=self.analyze_all).pack(side="right", padx=8)

        self.info = ttk.Label(self, textvariable=self.status_var, anchor="w")
        self.info.pack(fill="x", padx=12, pady=4)
        body = ttk.PanedWindow(self, orient="vertical")
        body.pack(fill="both", expand=True, padx=10, pady=5)
        table_frame = ttk.Frame(body)
        gel_frame = ttk.LabelFrame(body, text="Simulated electrophoresis gel (smaller fragments migrate farther)")
        body.add(table_frame, weight=3)
        body.add(gel_frame, weight=2)
        columns = ("combination", "sites", "cuts", "fragments")
        self.table = ttk.Treeview(table_frame, columns=columns, show="headings")
        headings = {"combination": "Enzyme combination", "sites": "Cut-site count", "cuts": "Cut positions (1-based)", "fragments": "Fragment lengths (bp)"}
        widths = {"combination": 210, "sites": 110, "cuts": 300, "fragments": 360}
        for col in columns:
            self.table.heading(col, text=headings[col])
            self.table.column(col, width=widths[col], anchor="w")
        scroll = ttk.Scrollbar(table_frame, orient="vertical", command=self.table.yview)
        self.table.configure(yscrollcommand=scroll.set)
        self.table.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")
        self.gel = tk.Canvas(gel_frame, background="#11151b", height=190, highlightthickness=0)
        self.gel.pack(fill="both", expand=True, padx=8, pady=8)
        self.gel.bind("<Configure>", lambda _event: self.draw_gel())

    def browse(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa *.fasta *.fna *.fas *.gz"), ("All files", "*.*")])
        if path:
            self.path_var.set(path)
            self.load()

    def load(self) -> None:
        try:
            self.description, self.sequence = read_fasta(self.path_var.get())
        except (OSError, ValueError) as exc:
            messagebox.showerror("Cannot load FASTA", str(exc))
            return
        self.status_var.set(f"Loaded {self.description}: {len(self.sequence):,} bp. Select enzymes and analyze.")
        self.analyze_all()

    def analyze_all(self) -> None:
        if not self.sequence:
            messagebox.showinfo("FASTA required", "Load an NCBI FASTA file first.")
            return
        selected = [name for name in ENZYMES if self.vars[name].get()]
        if not selected:
            messagebox.showinfo("Enzyme required", "Select at least one enzyme.")
            return
        self.table.delete(*self.table.get_children())
        for size in range(1, len(selected) + 1):
            for combination in itertools.combinations(selected, size):
                cuts, fragments = digest(self.sequence, combination)
                self.table.insert("", "end", values=(" + ".join(combination), len(cuts), ", ".join(map(str, cuts)) or "none", ", ".join(map(str, fragments))))
        self.status_var.set(f"Analyzed {sum(1 for size in range(1, len(selected)+1) for _ in itertools.combinations(selected, size))} combination(s).")
        self.draw_gel()

    def draw_gel(self) -> None:
        self.gel.delete("all")
        if not self.sequence:
            return
        selected = [name for name in ENZYMES if self.vars[name].get()]
        if not selected:
            return
        width, height = max(self.gel.winfo_width(), 500), max(self.gel.winfo_height(), 180)
        left, right, top, bottom = 90, 20, 20, 28
        self.gel.create_text(10, height - 12, text="larger fragments", fill="#cbd5e1", anchor="w")
        self.gel.create_text(width - 10, height - 12, text="smaller fragments", fill="#cbd5e1", anchor="e")
        combos = [tuple(c) for n in range(1, len(selected)+1) for c in itertools.combinations(selected, n)]
        for row, combination in enumerate(combos):
            _, fragments = digest(self.sequence, combination)
            y = top + (row + 0.5) * (height - top - bottom) / len(combos)
            self.gel.create_text(4, y, text="+".join(combination), fill="#f8fafc", anchor="w", font=("Segoe UI", 8))
            max_len = max(fragments)
            for frag in fragments:
                # Log scale gives useful separation for real genomes.
                x = left + (width - left - right) * (1 - (0.25 + 0.75 * (1 - min(frag / max_len, 1)) ** 0.55))
                self.gel.create_oval(x - 3, y - 3, x + 3, y + 3, fill="#65d9ff", outline="")


if __name__ == "__main__":
    DigestApp().mainloop()
