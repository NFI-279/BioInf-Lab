"""Lab 1, exercise 1: nucleotide composition."""
from dna_analysis import base_frequencies, clean_dna

EXAMPLE = 'ATGGCCAATGAGTGGAATCG'


def main() -> None:
    sequence = clean_dna(EXAMPLE)
    frequencies = base_frequencies(sequence)
    print('EXERCISE 1 - NUCLEOTIDE FREQUENCIES')
    print('DNA:', sequence)
    print('Sequence length:', len(sequence))
    print('\nBase   Count   Percentage')
    for base, count in frequencies.items():
        print(f'{base:<6} {count:<7} {100 * count / len(sequence):.2f}%')
    print('Total:', sum(frequencies.values()))


if __name__ == '__main__':
    main()
