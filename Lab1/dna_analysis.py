"""Shared, dependency-free DNA sequence analysis for Bioinformatics Lab 1."""
from collections import Counter

DNA_ALPHABET = frozenset('ACGT')


def clean_dna(sequence: str) -> str:
    """Normalize input and reject characters outside A, C, G, T."""
    normalized = ''.join(sequence.split()).upper()
    if not normalized:
        raise ValueError('DNA sequence cannot be empty')
    invalid = set(normalized) - DNA_ALPHABET
    if invalid:
        raise ValueError(f'Invalid DNA symbols: {", ".join(sorted(invalid))}')
    return normalized


def base_frequencies(sequence: str) -> dict[str, int]:
    """Count each nucleotide, including zero counts."""
    dna = clean_dna(sequence)
    counts = Counter(dna)
    return {base: counts[base] for base in 'ACGT'}


def count_kmers(sequence: str, size: int) -> dict[str, int]:
    """Count every overlapping DNA word of the requested length."""
    dna = clean_dna(sequence)
    if not 1 <= size <= len(dna):
        raise ValueError('Word size must be between 1 and the DNA length')
    counts = Counter(dna[position:position + size] for position in range(len(dna) - size + 1))
    return dict(sorted(counts.items()))
