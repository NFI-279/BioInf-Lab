"""Needleman-Wunsch global DNA alignment, implemented from first principles."""
from dataclasses import dataclass
from pathlib import Path
import gzip

DNA = set('ACGTRYSWKMBDHVN')


@dataclass
class AlignmentResult:
    sequence1: str
    sequence2: str
    scores: list[list[int]]
    directions: list[list[str]]
    path: list[tuple[int, int]]
    aligned1: str
    aligned2: str
    score: int
    matches: int
    mismatches: int
    gaps: int


def validate_dna(sequence: str) -> str:
    sequence = ''.join(sequence.split()).upper().replace('U', 'T')
    if not sequence:
        raise ValueError('DNA sequence cannot be empty.')
    invalid = sorted(set(sequence) - DNA)
    if invalid:
        raise ValueError('Invalid DNA character(s): ' + ', '.join(invalid))
    return sequence


def needleman_wunsch(seq1: str, seq2: str, match: int = 1,
                     mismatch: int = -1, gap: int = -2) -> AlignmentResult:
    """Global alignment; ties prefer diagonal, then up, then left."""
    seq1, seq2 = validate_dna(seq1), validate_dna(seq2)
    n, m = len(seq1), len(seq2)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    directions = [[''] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = i * gap
        directions[i][0] = '↑'
    for j in range(1, m + 1):
        scores[0][j] = j * gap
        directions[0][j] = '←'
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            up = scores[i - 1][j] + gap
            left = scores[i][j - 1] + gap
            best = max(diagonal, up, left)
            scores[i][j] = best
            directions[i][j] = '↖' if diagonal == best else ('↑' if up == best else '←')

    i, j = n, m
    path = [(i, j)]
    a1, a2 = [], []
    while i > 0 or j > 0:
        direction = directions[i][j]
        if direction == '↖':
            a1.append(seq1[i - 1]); a2.append(seq2[j - 1])
            i -= 1; j -= 1
        elif direction == '↑':
            a1.append(seq1[i - 1]); a2.append('-')
            i -= 1
        elif direction == '←':
            a1.append('-'); a2.append(seq2[j - 1])
            j -= 1
        else:
            raise RuntimeError(f'Invalid traceback at ({i}, {j})')
        path.append((i, j))
    a1, a2 = ''.join(reversed(a1)), ''.join(reversed(a2))
    matches = sum(x == y for x, y in zip(a1, a2) if x != '-' and y != '-')
    mismatches = sum(x != y for x, y in zip(a1, a2) if x != '-' and y != '-')
    gaps = sum(x == '-' or y == '-' for x, y in zip(a1, a2))
    return AlignmentResult(seq1, seq2, scores, directions, path, a1, a2,
                           scores[n][m], matches, mismatches, gaps)


def read_fasta(path: str) -> tuple[str, str, int]:
    """Read the first FASTA record, including .gz compressed files.

    Returns (header, sequence, number_of_records).
    """
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    header, chunks, count = None, [], 0
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('>'):
                count += 1
                if count == 1:
                    header = line[1:]
                continue
            if count == 1:
                chunks.append(line)
            elif count == 0:
                raise ValueError('Not a FASTA file: missing >header.')
    if header is None or not chunks:
        raise ValueError('No DNA record found in FASTA file.')
    return header, validate_dna(''.join(chunks)), count


def format_alignment(result: AlignmentResult, width: int = 70) -> str:
    a1, a2 = result.aligned1, result.aligned2
    lines = []
    for start in range(0, len(a1), width):
        x, y = a1[start:start + width], a2[start:start + width]
        marker = ''.join('|' if a == b and a != '-' else (' ' if '-' in (a, b) else '.')
                         for a, b in zip(x, y))
        lines.extend([f'S1  {x}', f'    {marker}', f'S2  {y}', ''])
    return '\n'.join(lines)
