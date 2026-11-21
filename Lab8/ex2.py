"""Lab 8 / Exercise 2: DNA sliding-window patterns and centers.

Run: python ex2.py
Requires matplotlib. Images are written into the outputs directory.
"""
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt

DEMO_DNA = 'CGGACTGATCTATCT' + ('A' * 27) + 'CGTAGCATCTATCGATCTATCTAGCGATCTA'
WINDOW = 30
OUTPUT = Path(__file__).resolve().parent / 'outputs'


def validate(sequence: str) -> str:
    clean = ''.join(sequence.upper().split())
    if not clean or set(clean) - set('ACGT'):
        raise ValueError('DNA sequence must consist only of A, C, G and T')
    return clean


def gc_percent(fragment: str) -> float:
    fragment = validate(fragment)
    return 100.0 * sum(letter in 'GC' for letter in fragment) / len(fragment)


def information_content(fragment: str) -> float:
    """Literal corrected whiteboard formula, NOT Shannon entropy.

    IC = 100 - mean over lags u=1..n-1 of
         (number of equal pairs at lag u) / ((n-u)*100).
    """
    fragment = validate(fragment)
    n = len(fragment)
    if n < 2:
        raise ValueError('Information content needs at least 2 bases')
    fractions = []
    for distance in range(1, n):
        equal_pairs = sum(left == right for left, right in
                          zip(fragment[:-distance], fragment[distance:]))
        fractions.append(equal_pairs / ((n - distance) * 100.0))
    return 100.0 - sum(fractions) / len(fractions)


@dataclass
class Pattern:
    name: str
    sequence: str
    gc: list[float]
    ic: list[float]
    center: tuple[float, float]


def calculate_pattern(name: str, sequence: str, window: int = WINDOW) -> Pattern:
    sequence = validate(sequence)
    if window < 2 or window > len(sequence):
        raise ValueError('Window must be between 2 and the sequence length')
    fragments = (sequence[i:i + window] for i in range(len(sequence) - window + 1))
    points = [(gc_percent(part), information_content(part)) for part in fragments]
    xs, ys = [x for x, _ in points], [y for _, y in points]
    # Arithmetic centroid of the plotted GC/IC points.
    center = (sum(xs) / len(xs), sum(ys) / len(ys))
    return Pattern(name, sequence, xs, ys, center)


def load_fasta(path: Path) -> str:
    import gzip
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if not lines or not lines[0].startswith('>'):
        raise ValueError('Not a FASTA file')
    body = []
    for line in lines[1:]:
        if line.startswith('>'): break  # first record only
        body.append(line)
    return validate(''.join(body))


def plot_pattern(pattern: Pattern, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(pattern.gc, pattern.ic, color='#6b50a1', linewidth=1.6, alpha=.85)
    ax.scatter(pattern.gc, pattern.ic, c=range(1, len(pattern.gc) + 1),
               cmap='viridis', s=22, zorder=3)
    ax.scatter(*pattern.center, marker='*', color='#d04b45', s=230,
               edgecolors='black', linewidths=.5, label='Pattern centroid', zorder=5)
    ax.set(xlabel='GC content (%)', ylabel='Information content',
           title=f'DNA sliding-window pattern | {pattern.name}')
    ax.grid(alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def plot_centers(patterns: list[Pattern], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 6))
    for pattern in patterns:
        x, y = pattern.center
        ax.scatter(x, y, marker='D', s=90)
        ax.annotate(pattern.name, (x, y), xytext=(8, 7), textcoords='offset points')
    ax.set(xlabel='Mean GC content (%)', ylabel='Mean information content',
           title='DNA promoter pattern centers')
    ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description='DNA promoter patterns (Lab 8)')
    parser.add_argument('--fasta', type=Path, action='append', default=[],
                        help='Optional FASTA promoter file (repeat for multiple sequences)')
    parser.add_argument('--window', type=int, default=WINDOW)
    args = parser.parse_args()
    patterns = [calculate_pattern('Example sequence', DEMO_DNA, args.window)]
    for path in args.fasta:
        patterns.append(calculate_pattern(path.stem, load_fasta(path), args.window))
    OUTPUT.mkdir(exist_ok=True)
    for pattern in patterns:
        print(f'{pattern.name}: {len(pattern.sequence)} bases, {len(pattern.gc)} windows')
        print(f'  Whole-sequence GC: {gc_percent(pattern.sequence):.2f}%')
        print(f'  Center: ({pattern.center[0]:.4f}, {pattern.center[1]:.4f})')
    plot_pattern(patterns[0], OUTPUT / 'dna_pattern.png')
    plot_centers(patterns, OUTPUT / 'promoter_centers.png')
    print('Saved images to:', OUTPUT)


if __name__ == '__main__':
    main()
