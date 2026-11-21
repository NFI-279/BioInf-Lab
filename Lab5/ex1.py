"""Lab 5, Exercise 1: generate 10 fragments and simulate gel electrophoresis."""
from pathlib import Path
import argparse
import math
import random
from dna_io import load_dna

ROOT = Path(__file__).resolve().parent
DEFAULT_FASTA = ROOT / 'data' / '480540219_genomic.fna'
MARKER_SIZES = [100, 200, 300, 400, 500, 700, 1000, 1500, 2000]

def fragments_from_dna(sequence, seed=42, count=10):
    region = sequence[:2000]
    if len(region) < 100:
        raise ValueError('The DNA input needs at least 100 bases')
    rng = random.Random(seed)
    fragments = []
    for _ in range(count):
        size = rng.randrange(100, len(region) + 1)
        offset = rng.randrange(len(region) - size + 1)
        fragments.append((offset + 1, region[offset:offset + size]))
    return fragments

def weight_da(base_pairs):
    return 617.96 * base_pairs + 36.04

def travel_distance(base_pairs, gel_length=10.0):
    if base_pairs <= 0:
        raise ValueError('Base-pair length must be positive')
    position = gel_length * math.log10(20000 / base_pairs) / math.log10(20000 / 100)
    return min(gel_length, max(0.0, position))

def draw_gel(samples, filename):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, ax = plt.subplots(figsize=(12, 7.5))
    ax.set_facecolor('#12213a')
    for lane in range(11):
        ax.add_patch(Rectangle((lane - .30, -.30), .60, .24, color='#687589'))
    for size in MARKER_SIZES:
        y = travel_distance(size)
        ax.plot([-.28, .28], [y, y], linewidth=5, color='#f8dd6b')
        ax.text(-.40, y, f'{size} bp', ha='right', va='center', color='white', fontsize=8)
    for lane, (_, fragment) in enumerate(samples, start=1):
        y = travel_distance(len(fragment))
        ax.plot([lane - .28, lane + .28], [y, y], linewidth=6, color='#6be5d0')
        ax.annotate(f'{len(fragment)} bp', (lane, y), xytext=(0, -15),
                    textcoords='offset points', ha='center', va='top', color='white', fontsize=8)
    ax.set_xlim(-1.4, 10.7)
    ax.set_ylim(10.6, -.6)
    ax.set_xticks(range(11), ['L'] + [str(i) for i in range(1, 11)])
    ax.set_xlabel('DNA samples (L = ladder)')
    ax.set_ylabel('Migration distance (cm)')
    ax.set_title('Lab 5 - Gel electrophoresis simulation (10 DNA fragments)')
    fig.tight_layout()
    fig.savefig(filename, dpi=180)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, default=DEFAULT_FASTA)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    title, dna = load_dna(args.fasta)
    samples = fragments_from_dna(dna, seed=args.seed)
    out = ROOT / 'outputs'
    out.mkdir(exist_ok=True)
    rows = ['Exercise 1 - Gel electrophoresis', f'Source: {args.fasta.name}',
            f'FASTA record: {title}', f'Full sequence length: {len(dna)} bp',
            'Approximate double-stranded DNA molecular weight:', '']
    for index, (start, fragment) in enumerate(samples, 1):
        bp = len(fragment)
        rows.append(f'Sample {index:02d}: start={start}, length={bp} bp, '
                    f'MW={weight_da(bp) / 1000:.2f} kDa, distance={travel_distance(bp):.2f} cm')
    (out / 'gel_results.txt').write_text('\n'.join(rows) + '\n', encoding='utf-8')
    draw_gel(samples, out / 'gel_electrophoresis.png')
    print('\n'.join(rows))
    print('\nSaved outputs/gel_electrophoresis.png and outputs/gel_results.txt')

if __name__ == '__main__':
    main()
