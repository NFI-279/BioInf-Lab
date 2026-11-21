"""Bioinformatics lab: simulate and detect artificial DNA transposons.

Simplified model (not a real-genome transposon caller):
- 200–400 bases TOTAL, including inserted elements.
- Each element has a 5-base terminal inverted repeat (TIR) at each end.
- A 4-base target-site duplication (TSD) directly flanks the element.
- Interior is 12–20 random bases, with no constraints on biological function.
- Detection scans for TSD + TIR + interior + reverse-complement(TIR) + TSD.
"""
import random
import tkinter as tk
from tkinter import ttk, messagebox

DNA = 'ACGT'
TIR_LENGTH = 5
TSD_LENGTH = 4
MIN_INTERIOR = 12
MAX_INTERIOR = 20
MIN_DNA = 200
MAX_DNA = 400


def random_dna(length, rng):
    return ''.join(rng.choice(DNA) for _ in range(length))


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def make_element(rng):
    tir = random_dna(TIR_LENGTH, rng)
    interior = random_dna(rng.randint(MIN_INTERIOR, MAX_INTERIOR), rng)
    return tir + interior + reverse_complement(tir)


def detect_transposons(sequence):
    """Find candidate element intervals from sequence alone; coordinates are 0-based.

    TSD flanks are outside the reported element interval. Because short DNA
    repeats occur by chance, a small number of false candidates is possible.
    """
    candidates = []
    n = len(sequence)
    for start in range(TSD_LENGTH, n - TSD_LENGTH):
        left_tsd = sequence[start - TSD_LENGTH:start]
        left_tir = sequence[start:start + TIR_LENGTH]
        if len(left_tir) != TIR_LENGTH:
            continue
        for interior_len in range(MIN_INTERIOR, MAX_INTERIOR + 1):
            end = start + 2 * TIR_LENGTH + interior_len
            if end + TSD_LENGTH > n:
                continue
            right_tir = sequence[end - TIR_LENGTH:end]
            right_tsd = sequence[end:end + TSD_LENGTH]
            if right_tir == reverse_complement(left_tir) and left_tsd == right_tsd:
                candidates.append({
                    'start': start, 'end': end,
                    'length': end - start,
                    'tir': left_tir, 'tsd': left_tsd,
                    'sequence': sequence[start:end],
                })
    # Prefer longer candidates when overlaps occur, then display in DNA order.
    selected = []
    for candidate in sorted(candidates, key=lambda x: (-x['length'], x['start'])):
        if all(candidate['end'] <= other['start'] or candidate['start'] >= other['end']
               for other in selected):
            selected.append(candidate)
    return sorted(selected, key=lambda x: x['start'])


def simulate(total_length, count, rng):
    """Create total-length sequence with nonoverlapping inserted elements.

    For reproducible, easy-to-read lab examples, random spacer lengths are >= 8.
    Each element has a duplicated 4-base target site on either side.
    """
    elements = [make_element(rng) for _ in range(count)]
    minimum_spacers = 8 * (count + 1)
    required = sum(len(element) + 2 * TSD_LENGTH for element in elements) + minimum_spacers
    if required > total_length:
        raise ValueError('The chosen DNA length is too short for this configuration.')

    extra = total_length - required
    spacer_lengths = [8] * (count + 1)
    for _ in range(extra):
        spacer_lengths[rng.randrange(count + 1)] += 1

    pieces = [random_dna(spacer_lengths[0], rng)]
    inserted = []
    position = len(pieces[0])
    for i, element in enumerate(elements):
        tsd = random_dna(TSD_LENGTH, rng)
        pieces.append(tsd + element + tsd)
        inserted.append({'start': position + TSD_LENGTH,
                         'end': position + TSD_LENGTH + len(element),
                         'sequence': element, 'tir': element[:TIR_LENGTH], 'tsd': tsd})
        position += len(pieces[-1])
        spacer = random_dna(spacer_lengths[i + 1], rng)
        pieces.append(spacer)
        position += len(spacer)
    sequence = ''.join(pieces)
    assert len(sequence) == total_length
    return sequence, inserted


class TransposonApp:
    def __init__(self, root):
        self.root = root
        self.root.title('Bioinformatics Lab — DNA Transposon Simulator')
        self.root.geometry('1050x760')
        self.root.minsize(800, 600)
        self.sequence = ''
        self.inserted = []
        self.detected = []

        top = ttk.Frame(root, padding=12)
        top.pack(fill='x')
        ttk.Label(top, text='DNA length (200–400):').pack(side='left')
        self.length_var = tk.StringVar(value='300')
        ttk.Entry(top, textvariable=self.length_var, width=7).pack(side='left', padx=(5, 18))
        ttk.Label(top, text='Transposons:').pack(side='left')
        self.count_var = tk.StringVar(value='4')
        ttk.Combobox(top, textvariable=self.count_var, values=('3', '4'), width=4,
                     state='readonly').pack(side='left', padx=(5, 18))
        ttk.Button(top, text='Generate DNA', command=self.generate).pack(side='left', padx=5)
        ttk.Button(top, text='Detect transposons', command=self.detect).pack(side='left', padx=5)

        ttk.Label(root, text='DNA sequence (highlighted regions = detected elements)',
                  padding=(12, 3)).pack(anchor='w')
        sequence_frame = ttk.Frame(root, padding=(12, 0, 12, 10))
        sequence_frame.pack(fill='both', expand=True)
        self.sequence_text = tk.Text(sequence_frame, height=13, wrap='char',
                                     font=('Consolas', 12), spacing2=3)
        self.sequence_text.pack(side='left', fill='both', expand=True)
        scrollbar = ttk.Scrollbar(sequence_frame, orient='vertical',
                                  command=self.sequence_text.yview)
        scrollbar.pack(side='right', fill='y')
        self.sequence_text.configure(yscrollcommand=scrollbar.set)
        self.sequence_text.tag_configure('element', background='#ffe69a', foreground='#242424')
        self.sequence_text.configure(state='disabled')

        ttk.Label(root, text='Detected elements (1-based positions, inclusive)',
                  padding=(12, 3)).pack(anchor='w')
        table_frame = ttk.Frame(root, padding=(12, 0, 12, 8))
        table_frame.pack(fill='both', expand=True)
        columns = ('id', 'start', 'end', 'length', 'tir', 'tsd', 'verified')
        self.table = ttk.Treeview(table_frame, columns=columns, show='headings', height=8)
        headers = ('#', 'Start', 'End', 'Length', 'Inverted repeat', 'Direct repeat', 'Matches insertion')
        widths = (45, 75, 75, 75, 145, 130, 140)
        for col, header, width in zip(columns, headers, widths):
            self.table.heading(col, text=header)
            self.table.column(col, width=width, anchor='center')
        self.table.pack(side='left', fill='both', expand=True)
        table_scroll = ttk.Scrollbar(table_frame, orient='vertical', command=self.table.yview)
        table_scroll.pack(side='right', fill='y')
        self.table.configure(yscrollcommand=table_scroll.set)
        self.table.bind('<<TreeviewSelect>>', self.select_element)

        self.status_var = tk.StringVar(value='Generate a DNA sequence to begin.')
        ttk.Label(root, textvariable=self.status_var, padding=12, wraplength=980).pack(fill='x')
        ttk.Label(root, text='Model: 5-base terminal inverted repeats + 4-base flanking direct repeats (target-site duplication).',
                  padding=(12, 0, 12, 12)).pack(anchor='w')

    def set_sequence(self, sequence):
        self.sequence_text.configure(state='normal')
        self.sequence_text.delete('1.0', 'end')
        self.sequence_text.insert('1.0', sequence)
        self.sequence_text.tag_remove('element', '1.0', 'end')
        self.sequence_text.configure(state='disabled')

    def generate(self):
        try:
            length = int(self.length_var.get())
            count = int(self.count_var.get())
            if not MIN_DNA <= length <= MAX_DNA:
                raise ValueError('DNA length must be between 200 and 400 bases.')
            if count not in (3, 4):
                raise ValueError('Choose 3 or 4 transposons.')
            rng = random.Random()
            for _ in range(1000):
                sequence, inserted = simulate(length, count, rng)
                expected = {(e['start'], e['end']) for e in inserted}
                actual = {(e['start'], e['end']) for e in detect_transposons(sequence)}
                if actual == expected:
                    break
            else:
                raise ValueError('Could not generate an unambiguous example; please retry.')
            self.sequence, self.inserted = sequence, inserted
        except ValueError as exc:
            messagebox.showerror('Invalid input', str(exc))
            return
        self.detected = []
        self.set_sequence(self.sequence)
        self.table.delete(*self.table.get_children())
        self.status_var.set(f'Generated {len(self.sequence)} bases with {count} inserted elements. Click Detect transposons.')

    def detect(self):
        if not self.sequence:
            messagebox.showinfo('No DNA', 'Generate a DNA sequence first.')
            return
        self.detected = detect_transposons(self.sequence)
        expected = {(e['start'], e['end']) for e in self.inserted}
        actual = {(e['start'], e['end']) for e in self.detected}
        self.table.delete(*self.table.get_children())
        self.sequence_text.configure(state='normal')
        self.sequence_text.tag_remove('element', '1.0', 'end')
        for i, element in enumerate(self.detected, start=1):
            start, end = element['start'], element['end']
            self.table.insert('', 'end', iid=str(i - 1), values=(
                i, start + 1, end, element['length'], element['tir'],
                element['tsd'], 'Yes' if (start, end) in expected else 'No'))
            self.sequence_text.tag_add('element', f'1.0+{start}c', f'1.0+{end}c')
        self.sequence_text.configure(state='disabled')
        matched = len(expected & actual)
        self.status_var.set(f'Detected {len(actual)} candidate(s); {matched}/{len(expected)} inserted elements found. '
                            f'False positives: {len(actual - expected)}. '
                            'Coordinates are 1-based and inclusive.')

    def select_element(self, _event):
        selection = self.table.selection()
        if not selection:
            return
        element = self.detected[int(selection[0])]
        self.sequence_text.see(f"1.0+{element['start']}c")


def main():
    root = tk.Tk()
    TransposonApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
