"""Lab 8 / Exercise 1: global alignment with Needleman–Wunsch.

Run: python ex1.py
Requires matplotlib. Images are written into the outputs directory.
"""
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

FIRST = 'ACCGTGAAGCCAATAC'
SECOND = 'AGCGTGCAGCCAATAC'
OUTPUT = Path(__file__).resolve().parent / 'outputs'


@dataclass
class Alignment:
    first: str
    second: str
    score: int
    grid: list[list[int]]
    route: list[tuple[int, int]]


def align(a: str, b: str, match: int = 1, mismatch: int = -1, gap: int = -2) -> Alignment:
    """Build the DP table and follow the chosen predecessors back to (0, 0)."""
    a, b = a.upper(), b.upper()
    if not a or not b or any(ch not in 'ACGT' for ch in a + b):
        raise ValueError('Both sequences must contain only A, C, G, T and be nonempty')
    height, width = len(a) + 1, len(b) + 1
    grid = [[0] * width for _ in range(height)]
    parent = [[''] * width for _ in range(height)]
    for r in range(1, height):
        grid[r][0], parent[r][0] = r * gap, 'U'
    for c in range(1, width):
        grid[0][c], parent[0][c] = c * gap, 'L'
    for r, letter_a in enumerate(a, start=1):
        for c, letter_b in enumerate(b, start=1):
            candidates = (
                (grid[r - 1][c - 1] + (match if letter_a == letter_b else mismatch), 'D'),
                (grid[r - 1][c] + gap, 'U'),
                (grid[r][c - 1] + gap, 'L'),
            )
            grid[r][c], parent[r][c] = max(candidates, key=lambda item: item[0])
    r, c = len(a), len(b)
    backwards_a, backwards_b = [], []
    route = [(r, c)]
    while r or c:
        step = parent[r][c]
        if step == 'D':
            backwards_a.append(a[r - 1]); backwards_b.append(b[c - 1]); r -= 1; c -= 1
        elif step == 'U':
            backwards_a.append(a[r - 1]); backwards_b.append('-'); r -= 1
        else:
            backwards_a.append('-'); backwards_b.append(b[c - 1]); c -= 1
        route.append((r, c))
    return Alignment(''.join(reversed(backwards_a)), ''.join(reversed(backwards_b)), grid[-1][-1], grid, route)


def statistics(result: Alignment) -> dict[str, int]:
    pairs = zip(result.first, result.second)
    count = {'matches': 0, 'mismatches': 0, 'gaps': 0}
    for x, y in pairs:
        if '-' in (x, y): count['gaps'] += 1
        elif x == y: count['matches'] += 1
        else: count['mismatches'] += 1
    return count


def draw_scoring_matrix(a: str, b: str, result: Alignment, destination: Path) -> None:
    size = max(len(a), len(b))
    fig, ax = plt.subplots(figsize=(max(8, size * .65), max(7, size * .58)))
    values = result.grid
    low, high = min(map(min, values)), max(map(max, values))
    img = ax.imshow(values, cmap='PuOr', vmin=low, vmax=high)
    ax.set_xticks(range(len(b) + 1), ['–', *b])
    ax.set_yticks(range(len(a) + 1), ['–', *a])
    ax.set_xticks([i - .5 for i in range(len(b) + 2)], minor=True)
    ax.set_yticks([i - .5 for i in range(len(a) + 2)], minor=True)
    ax.grid(which='minor', color='white', linewidth=.7)
    ax.tick_params(which='minor', bottom=False, left=False)
    for row, cells in enumerate(values):
        for col, value in enumerate(cells):
            ax.text(col, row, str(value), ha='center', va='center', fontsize=8,
                    color='white' if value < (low + high) / 2 else 'black')
    path_x = [col for row, col in reversed(result.route)]
    path_y = [row for row, col in reversed(result.route)]
    ax.plot(path_x, path_y, '-o', color='#00dfaa', linewidth=2.2, markersize=4,
            label='Optimal traceback')
    ax.set(title='Needleman–Wunsch | Score matrix and traceback', xlabel='Sequence 2', ylabel='Sequence 1')
    ax.legend(loc='upper right')
    fig.colorbar(img, ax=ax, shrink=.8, label='DP score')
    fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)


def draw_alignment(result: Alignment, destination: Path) -> None:
    n = len(result.first)
    fig, ax = plt.subplots(figsize=(max(10, n * .65), 3.5))
    ax.set_xlim(-.7, n - .3); ax.set_ylim(-.5, 2.7)
    ax.set_yticks([2, 0], ['S1', 'S2']); ax.set_xticks(range(n), range(1, n + 1))
    ax.tick_params(length=0)
    for i, (x, y) in enumerate(zip(result.first, result.second)):
        tone = '#d9f2e6' if x == y else ('#ffe9c7' if '-' in (x, y) else '#f7d9dd')
        ax.axvspan(i - .48, i + .48, ymin=.18, ymax=.85, color=tone)
        ax.text(i, 2, x, ha='center', va='center', family='monospace', weight='bold')
        ax.text(i, 1, '|' if x == y else (' ' if '-' in (x, y) else '×'), ha='center', va='center')
        ax.text(i, 0, y, ha='center', va='center', family='monospace', weight='bold')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.set(title=f'Global DNA alignment | Optimal score: {result.score}', xlabel='Alignment position')
    fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    result = align(FIRST, SECOND)
    symbols = ''.join('|' if x == y else (' ' if '-' in (x, y) else '.')
                      for x, y in zip(result.first, result.second))
    print(result.first, symbols, result.second, sep='\n')
    print('Score:', result.score, '|', statistics(result))
    draw_scoring_matrix(FIRST, SECOND, result, OUTPUT / 'scoring_matrix.png')
    draw_alignment(result, OUTPUT / 'dna_alignment.png')
    print('Saved images to:', OUTPUT)


if __name__ == '__main__':
    main()
