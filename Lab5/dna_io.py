"""Read the first DNA record in an ordinary or gzipped FASTA file."""
from pathlib import Path
import gzip

def load_dna(filename):
    path = Path(filename)
    open_file = gzip.open if path.name.lower().endswith('.gz') else open
    header = None
    fragments = []
    with open_file(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    break  # first FASTA record only
                header = line[1:]
            elif header is None:
                raise ValueError('Expected a FASTA header starting with >')
            else:
                fragments.append(line.upper())
    sequence = ''.join(fragments)
    if not sequence:
        raise ValueError('The FASTA record has no sequence')
    if set(sequence) - set('ACGTN'):
        raise ValueError('Unexpected character in DNA sequence')
    return header, sequence
