"""Exercise 2 - simulated transposomes made from inverted repeats.

An inverted-repeat transposome has this form::

    LEFT_REPEAT + spacer (at least 10 bases) + reverse_complement(LEFT_REPEAT)

Coordinates printed by this program are 1-based and inclusive.  The reported
length is the complete interval from the first base of the left repeat to the
last base of the right repeat.
"""
from __future__ import annotations

import argparse
import random
from dataclasses import dataclass

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Transposome:
    name: str
    start: int       # zero-based internally
    end: int         # exclusive internally
    repeat_length: int
    spacer_length: int

    @property
    def position(self) -> tuple[int, int]:
        return self.start + 1, self.end

    @property
    def length(self) -> int:
        return self.end - self.start


def interval_relation(first: Transposome, second: Transposome) -> str:
    """Classify two intervals, assuming they are different transposomes."""
    if first.end <= second.start or second.end <= first.start:
        return "normal/disjoint"
    contains = (first.start <= second.start and first.end >= second.end) or (
        second.start <= first.start and second.end >= first.end
    )
    return "nested" if contains else "overlapping"


def random_dna(length: int, rng: random.Random) -> list[str]:
    return [rng.choice(BASES) for _ in range(length)]


def insert_transposome(
    genome: list[str], start: int, repeat: str, spacer_length: int, name: str
) -> Transposome:
    """Insert an inverted-repeat pair, overwriting the generated background."""
    right = reverse_complement(repeat)
    end = start + len(repeat) + spacer_length + len(right)
    if start < 0 or end > len(genome):
        raise ValueError(f"{name} does not fit in the generated DNA sequence")
    genome[start : start + len(repeat)] = repeat
    genome[start + len(repeat) + spacer_length : end] = right
    return Transposome(name, start, end, len(repeat), spacer_length)


def build_demo(length: int = 260, seed: int = 42) -> tuple[str, list[Transposome]]:
    """Create examples of all three requested positional cases."""
    rng = random.Random(seed)
    genome = random_dna(length, rng)
    # The positions deliberately produce: disjoint, overlapping, and nested.
    specs = [
        ("T1", 10, "ACGTAC", 10),   # [11, 32], disjoint from T2
        ("T2", 45, "GATTACA", 10),  # [46, 69], normal/disjoint
        ("T3", 85, "CCGTA", 10),    # [86, 105], overlaps T4
        ("T4", 98, "TGCATG", 10),   # [99, 125], overlapping T3
        ("T5", 155, "AGTCGATA", 10),  # [156, 181], contains T6
        ("T6", 160, "CGAT", 10),    # [161, 188], nested in T5 (overwritten arms are intentional)
    ]
    transposomes = [insert_transposome(genome, start, repeat, gap, name) for name, start, repeat, gap in specs]
    return "".join(genome), transposomes


def print_report(sequence: str, transposomes: list[Transposome], seed: int) -> None:
    print(f"Random DNA length: {len(sequence)} bp | seed: {seed}")
    print("\nTransposome intervals (1-based, inclusive):")
    print("Name   Position       Length   Repeat  Spacer")
    for item in transposomes:
        start, end = item.position
        print(f"{item.name:<6} {start:>4}-{end:<7} {item.length:>6} {item.repeat_length:>7} {item.spacer_length:>7}")
    print("\nRelations:")
    for i, first in enumerate(transposomes):
        for second in transposomes[i + 1 :]:
            relation = interval_relation(first, second)
            if relation != "normal/disjoint" or {first.name, second.name} <= {"T1", "T2"}:
                print(f"{first.name} / {second.name}: {relation}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate and analyze inverted-repeat transposomes.")
    parser.add_argument("--length", type=int, default=260, help="Generated DNA length (default: 260)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible DNA")
    args = parser.parse_args()
    if args.length < 220:
        parser.error("--length must be at least 220 for the demonstration intervals")
    sequence, transposomes = build_demo(args.length, args.seed)
    print_report(sequence, transposomes, args.seed)


if __name__ == "__main__":
    main()
