"""Tkinter teaching GUI for Needleman-Wunsch DNA global alignment."""
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path

from needleman_wunsch import needleman_wunsch, read_fasta, format_alignment

DEFAULT_1 = 'ACCGTGAAGCCAATAC'
DEFAULT_2 = 'AGCGTGCAGCCAATAC'
MAX_LENGTH = 120  # 121x121 cells; deliberate cap for a readable GUI
CELL_W, CELL_H = 60, 33


class AlignmentApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Needleman–Wunsch | DNA Global Alignment Lab')
        self.geometry('1200x850')
        self.minsize(850, 620)
        self.configure(bg='#f2f5f9')
        self.result = None
        self.trace_index = 0
        self.source_sequences = [DEFAULT_1, DEFAULT_2]
        self.source_names = ['Manual sequence', 'Manual sequence']
        self.entries = []
        self.start_vars = []
        self.length_vars = []
        self.info_vars = []
        self._build()
        self.align()

    def _build(self):
        style = ttk.Style(self)
        if 'clam' in style.theme_names():
            style.theme_use('clam')
        style.configure('TNotebook.Tab', padding=(14, 7))
        style.configure('Action.TButton', padding=(10, 7))
        top = ttk.Frame(self, padding=12)
        top.pack(fill='x')
        ttk.Label(top, text='Needleman–Wunsch DNA Alignment', font=('Arial', 18, 'bold')).pack(anchor='w')
        ttk.Label(top, text='Global alignment • dynamic programming • complete scoring matrix and traceback').pack(anchor='w', pady=(3, 10))

        sequences = ttk.LabelFrame(top, text='DNA sequences / FASTA input', padding=8)
        sequences.pack(fill='x')
        for index, default in enumerate((DEFAULT_1, DEFAULT_2)):
            row = ttk.Frame(sequences)
            row.pack(fill='x', pady=3)
            ttk.Label(row, text=f'S{index + 1}', width=4).pack(side='left')
            entry = ttk.Entry(row, width=58)
            entry.insert(0, default)
            entry.pack(side='left', fill='x', expand=True, padx=(0, 7))
            self.entries.append(entry)
            ttk.Button(row, text='Load FASTA…', command=lambda idx=index: self.load_fasta(idx)).pack(side='left')
            ttk.Label(row, text='Start:').pack(side='left', padx=(10, 3))
            start = tk.StringVar(value='1')
            self.start_vars.append(start)
            ttk.Entry(row, textvariable=start, width=5).pack(side='left')
            ttk.Label(row, text='Length:').pack(side='left', padx=(7, 3))
            length = tk.StringVar(value=str(MAX_LENGTH))
            self.length_vars.append(length)
            ttk.Entry(row, textvariable=length, width=5).pack(side='left')
            info = tk.StringVar(value='Manual sequence (16 bases)')
            self.info_vars.append(info)
            ttk.Label(sequences, textvariable=info, foreground='#52657c').pack(anchor='w', padx=32)
        ttk.Label(sequences, text='Start and Length are used when loading FASTA. The text boxes contain the selected region and can be edited manually.',
                  foreground='#52657c').pack(anchor='w', pady=(5, 0))

        settings = ttk.Frame(top)
        settings.pack(fill='x', pady=(10, 0))
        self.match_var = tk.StringVar(value='1')
        self.mismatch_var = tk.StringVar(value='-1')
        self.gap_var = tk.StringVar(value='-2')
        for title, var in (('Match', self.match_var), ('Mismatch', self.mismatch_var), ('Gap', self.gap_var)):
            ttk.Label(settings, text=f'{title}:').pack(side='left', padx=(0, 4))
            ttk.Entry(settings, textvariable=var, width=6).pack(side='left', padx=(0, 14))
        ttk.Button(settings, text='Align sequences', style='Action.TButton', command=self.align).pack(side='left')
        ttk.Button(settings, text='Reset examples', command=self.reset).pack(side='left', padx=8)
        ttk.Button(settings, text='Save results…', command=self.save_results).pack(side='left')
        self.summary = tk.StringVar(value='')
        ttk.Label(top, textvariable=self.summary, font=('Arial', 11, 'bold'), foreground='#153d68').pack(anchor='w', pady=(10, 0))

        self.tabs = ttk.Notebook(self)
        self.tabs.pack(fill='both', expand=True, padx=12, pady=(0, 12))
        matrix_tab = ttk.Frame(self.tabs)
        traceback_tab = ttk.Frame(self.tabs, padding=10)
        alignment_tab = ttk.Frame(self.tabs, padding=10)
        self.tabs.add(matrix_tab, text='Scoring matrix')
        self.tabs.add(traceback_tab, text='Traceback steps')
        self.tabs.add(alignment_tab, text='Final alignment & algorithm')

        ttk.Label(matrix_tab, text='Rows = S1; columns = S2. Each cell shows score and the chosen direction. Gold = traceback; green = current step.',
                  padding=8).pack(anchor='w')
        frame = ttk.Frame(matrix_tab)
        frame.pack(fill='both', expand=True)
        self.canvas = tk.Canvas(frame, bg='white', highlightthickness=0)
        hscroll = ttk.Scrollbar(frame, orient='horizontal', command=self.canvas.xview)
        vscroll = ttk.Scrollbar(frame, orient='vertical', command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hscroll.set, yscrollcommand=vscroll.set)
        self.canvas.grid(row=0, column=0, sticky='nsew')
        vscroll.grid(row=0, column=1, sticky='ns')
        hscroll.grid(row=1, column=0, sticky='ew')
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        self.canvas.bind('<Control-MouseWheel>', lambda e: self.canvas.xview_scroll(-int(e.delta / 120), 'units'))
        self.canvas.bind('<MouseWheel>', lambda e: self.canvas.yview_scroll(-int(e.delta / 120), 'units'))
        controls = ttk.Frame(matrix_tab, padding=8)
        controls.pack(fill='x')
        ttk.Button(controls, text='◀ Previous traceback step', command=lambda: self.step_trace(-1)).pack(side='left')
        ttk.Button(controls, text='Next traceback step ▶', command=lambda: self.step_trace(1)).pack(side='left', padx=8)
        self.trace_label = tk.StringVar()
        ttk.Label(controls, textvariable=self.trace_label).pack(side='left', padx=10)

        ttk.Label(traceback_tab, text='Traceback starts in the bottom-right cell and ends at (0, 0).',
                  padding=(0, 0, 0, 8)).pack(anchor='w')
        columns = ('step', 'cell', 'move', 's1', 's2', 'meaning')
        self.trace_tree = ttk.Treeview(traceback_tab, columns=columns, show='headings', height=14)
        headings = ('Step', 'Cell (row, col)', 'Direction', 'S1 character', 'S2 character', 'Operation')
        for key, heading, width in zip(columns, headings, (65, 125, 100, 115, 115, 280)):
            self.trace_tree.heading(key, text=heading)
            self.trace_tree.column(key, width=width, anchor='center' if key != 'meaning' else 'w')
        self.trace_tree.pack(fill='both', expand=True)
        self.trace_tree.bind('<<TreeviewSelect>>', self.on_trace_selected)

        self.alignment_text = tk.Text(alignment_tab, font=('Consolas', 11), wrap='none', height=14)
        self.alignment_text.pack(fill='both', expand=True)
        ttk.Label(alignment_tab, text='Recurrence: F[i,j] = max(F[i−1,j−1] + substitution, F[i−1,j] + gap, F[i,j−1] + gap).\n'
                  'Initialize F[i,0] = i × gap and F[0,j] = j × gap. Ties: diagonal → up → left.\n'
                  'Traceback: ↖ aligns two letters, ↑ aligns S1 with a gap, ← aligns a gap with S2.',
                  justify='left', padding=(0, 12, 0, 0)).pack(anchor='w')

    def reset(self):
        for idx, seq in enumerate((DEFAULT_1, DEFAULT_2)):
            self.source_sequences[idx] = seq
            self.source_names[idx] = 'Manual sequence'
            self.start_vars[idx].set('1')
            self.length_vars[idx].set(str(MAX_LENGTH))
            self.entries[idx].delete(0, 'end')
            self.entries[idx].insert(0, seq)
            self.info_vars[idx].set('Manual sequence (16 bases)')
        self.match_var.set('1'); self.mismatch_var.set('-1'); self.gap_var.set('-2')
        self.align()

    def load_fasta(self, idx):
        filename = filedialog.askopenfilename(title=f'Load FASTA for S{idx + 1}',
                    filetypes=[('FASTA files', '*.fasta *.fa *.fna *.fas *.gz'), ('All files', '*.*')])
        if not filename:
            return
        try:
            header, sequence, records = read_fasta(filename)
            self.source_sequences[idx] = sequence
            self.source_names[idx] = f'{Path(filename).name}: {header}'
            self._apply_region(idx)
            suffix = f' (first of {records} records)' if records > 1 else ''
            self.info_vars[idx].set(f'{Path(filename).name}: {len(sequence):,} bases{suffix} | {header[:85]}')
        except (OSError, ValueError, EOFError) as exc:
            messagebox.showerror('FASTA loading error', str(exc))

    def _apply_region(self, idx):
        start = int(self.start_vars[idx].get())
        length = int(self.length_vars[idx].get())
        source = self.source_sequences[idx]
        if start < 1 or start > len(source):
            raise ValueError(f'S{idx + 1}: Start must be between 1 and {len(source)}.')
        if not 1 <= length <= MAX_LENGTH:
            raise ValueError(f'S{idx + 1}: Length must be 1–{MAX_LENGTH}.')
        selected = source[start - 1:start - 1 + length]
        self.entries[idx].delete(0, 'end')
        self.entries[idx].insert(0, selected)

    def align(self):
        try:
            seq1 = self.entries[0].get().strip()
            seq2 = self.entries[1].get().strip()
            if max(len(seq1), len(seq2)) > MAX_LENGTH:
                raise ValueError(f'This visualization supports up to {MAX_LENGTH} bases per sequence. Select a shorter FASTA region.')
            match = int(self.match_var.get())
            mismatch = int(self.mismatch_var.get())
            gap = int(self.gap_var.get())
            self.result = needleman_wunsch(seq1, seq2, match, mismatch, gap)
        except ValueError as exc:
            messagebox.showerror('Alignment error', str(exc))
            return
        r = self.result
        self.summary.set(f'Optimal global score: {r.score}    |    Matches: {r.matches}    |    Mismatches: {r.mismatches}    |    Gap columns: {r.gaps}    |    Matrix: {len(r.sequence1)+1} × {len(r.sequence2)+1}')
        self.trace_index = 0
        self.draw_matrix()
        self.show_trace_steps()
        self.show_alignment()

    def draw_matrix(self):
        r = self.result
        self.canvas.delete('all')
        rows, cols = len(r.scores), len(r.scores[0])
        path_cells = set(r.path)
        current = r.path[self.trace_index]
        # Header row and column occupy index 0 in canvas; DP cell (i,j) occupies (i+1,j+1).
        for i in range(rows + 1):
            for j in range(cols + 1):
                x0, y0 = j * CELL_W, i * CELL_H
                x1, y1 = x0 + CELL_W, y0 + CELL_H
                if i == 0 and j == 0:
                    label, bg = 'S1 / S2', '#dfe7f0'
                elif i == 0:
                    label, bg = ('–' if j == 1 else r.sequence2[j - 2]), '#dfe7f0'
                elif j == 0:
                    label, bg = ('–' if i == 1 else r.sequence1[i - 2]), '#dfe7f0'
                else:
                    pi, pj = i - 1, j - 1
                    direction = r.directions[pi][pj]
                    label = f'{r.scores[pi][pj]} {direction}'
                    bg = '#b8efcf' if (pi, pj) == current else ('#ffe8aa' if (pi, pj) in path_cells else ('#ffffff' if (pi + pj) % 2 == 0 else '#f5f8fc'))
                self.canvas.create_rectangle(x0, y0, x1, y1, fill=bg, outline='#d5dde7')
                self.canvas.create_text((x0 + x1) / 2, (y0 + y1) / 2, text=label,
                                        font=('Consolas', 9, 'bold' if i == 0 or j == 0 else 'normal'), fill='#243b53')
        self.canvas.configure(scrollregion=(0, 0, (cols + 1) * CELL_W, (rows + 1) * CELL_H))
        self.trace_label.set(f'Traceback step {self.trace_index + 1} / {len(r.path)}  |  cell {current}')

    def step_trace(self, delta):
        if not self.result:
            return
        self.trace_index = max(0, min(len(self.result.path) - 1, self.trace_index + delta))
        self.draw_matrix()
        i, j = self.result.path[self.trace_index]
        self.canvas.xview_moveto(max(0, (j * CELL_W - 200) / max(1, (len(self.result.sequence2) + 2) * CELL_W)))
        self.canvas.yview_moveto(max(0, (i * CELL_H - 100) / max(1, (len(self.result.sequence1) + 2) * CELL_H)))
        item = str(self.trace_index)
        if self.trace_tree.exists(item):
            self.trace_tree.selection_set(item)
            self.trace_tree.see(item)

    def show_trace_steps(self):
        for item in self.trace_tree.get_children():
            self.trace_tree.delete(item)
        r = self.result
        for index, (i, j) in enumerate(r.path):
            direction = r.directions[i][j]
            if direction == '↖':
                c1, c2 = r.sequence1[i - 1], r.sequence2[j - 1]
                meaning = 'Match' if c1 == c2 else 'Mismatch'
            elif direction == '↑':
                c1, c2, meaning = r.sequence1[i - 1], '–', 'Gap in S2'
            elif direction == '←':
                c1, c2, meaning = '–', r.sequence2[j - 1], 'Gap in S1'
            else:
                c1, c2, meaning = '–', '–', 'Finished at origin'
            self.trace_tree.insert('', 'end', iid=str(index), values=(index + 1, f'({i}, {j})', direction or 'START', c1, c2, meaning))

    def on_trace_selected(self, _event):
        selected = self.trace_tree.selection()
        if selected and self.result:
            index = int(selected[0])
            if index != self.trace_index:
                self.trace_index = index
                self.draw_matrix()

    def show_alignment(self):
        r = self.result
        text = (f'S1: {r.sequence1}\nS2: {r.sequence2}\n\n'
                f'Global alignment score: {r.score}\n'
                f'Matches: {r.matches} | Mismatches: {r.mismatches} | Gaps: {r.gaps}\n\n'
                + format_alignment(r))
        self.alignment_text.configure(state='normal')
        self.alignment_text.delete('1.0', 'end')
        self.alignment_text.insert('1.0', text)
        self.alignment_text.configure(state='disabled')

    def save_results(self):
        if not self.result:
            return
        path = filedialog.asksaveasfilename(title='Save alignment report', defaultextension='.txt',
                                            filetypes=[('Text report', '*.txt')])
        if not path:
            return
        try:
            r = self.result
            lines = ['Needleman–Wunsch DNA global alignment',
                     f'Scoring: match={self.match_var.get()}, mismatch={self.mismatch_var.get()}, gap={self.gap_var.get()}',
                     f'Score={r.score}, matches={r.matches}, mismatches={r.mismatches}, gaps={r.gaps}',
                     '', format_alignment(r), 'Score matrix (rows=S1, columns=S2):']
            for row in r.scores:
                lines.append('\t'.join(map(str, row)))
            lines.append('\nTraceback (from end to origin):')
            for i, j in r.path:
                lines.append(f'({i},{j}) {r.directions[i][j]}')
            Path(path).write_text('\n'.join(lines) + '\n', encoding='utf-8')
            messagebox.showinfo('Saved', f'Report saved to:\n{path}')
        except OSError as exc:
            messagebox.showerror('Save error', str(exc))


if __name__ == '__main__':
    AlignmentApp().mainloop()
