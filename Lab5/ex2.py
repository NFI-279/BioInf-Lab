"""Lab 5, Exercise 2: find overlapping DNA motifs of lengths 3, 4 and 5."""
from pathlib import Path
from collections import Counter
import argparse
from dna_io import load_dna

ROOT = Path(__file__).resolve().parent
DEFAULT_FASTA = ROOT / 'data' / 'GCF_000919155.1_ViralProj244477_genomic.fna.gz'

def motif_frequencies(dna, sizes=(3, 4, 5)):
    counts = Counter()
    for length in sizes:
        for start in range(len(dna) - length + 1):
            word = dna[start:start + length]
            if all(base in 'ACGT' for base in word):
                counts[word] += 1
    return counts

def matching_motifs(counts, minimum=6, maximum=10):
    if minimum > maximum:
        raise ValueError('Minimum cannot exceed maximum')
    return sorted(((word, count) for word, count in counts.items()
                   if minimum <= count <= maximum), key=lambda pair: (len(pair[0]), pair[0]))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', type=Path, default=DEFAULT_FASTA)
    parser.add_argument('--min', type=int, default=6)
    parser.add_argument('--max', type=int, default=10)
    args = parser.parse_args()
    header, dna = load_dna(args.fasta)
    found = matching_motifs(motif_frequencies(dna), args.min, args.max)
    output = ROOT / 'outputs'
    output.mkdir(exist_ok=True)
    lines = [f'{count} {word}' for word, count in found]
    (output / 'motif_frequencies.txt').write_text('\n'.join(lines) + ('\n' if lines else ''), encoding='utf-8')
    print(f'FASTA: {args.fasta.name} | {header}')
    print(f'Sequence length: {len(dna)} bases')
    print(f'3-5-base motifs appearing {args.min} to {args.max} times: {len(found)}')
    print('First results: ' + (', '.join(f'{word}={count}' for word, count in found[:12]) or '(none)'))
    print('Saved outputs/motif_frequencies.txt')

if __name__ == '__main__':
    main()
