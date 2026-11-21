"""Lab 1, exercise 2: overlapping dinucleotide and trinucleotide counts."""
from dna_analysis import clean_dna, count_kmers

EXAMPLE = 'ATGGCCAATGAGTGGAATCG' + 'C' * 20


def display_words(sequence: str, word_length: int) -> None:
    frequencies = count_kmers(sequence, word_length)
    print(f'\n{word_length}-NUCLEOTIDE COMBINATIONS')
    print('Word   Occurrences')
    for word, occurrences in frequencies.items():
        print(f'{word:<6} {occurrences}')
    print('Distinct combinations:', len(frequencies))
    print('Total overlapping windows:', sum(frequencies.values()))


def main() -> None:
    dna = clean_dna(EXAMPLE)
    print('EXERCISE 2 - DNA WORD FREQUENCIES')
    print('DNA:', dna)
    print('Sequence length:', len(dna))
    for size in (2, 3):
        display_words(dna, size)


if __name__ == '__main__':
    main()
