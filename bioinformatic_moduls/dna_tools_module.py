from typing import Tuple


def is_dna_or_rna(seq: str) -> bool:
    return 't' in seq.lower()


COMPLEMENT_DNA = {
    'A': 'T',
    'a': 't',
    'G': 'C',
    'g': 'c',

    'T': 'A',
    't': 'a',
    'C': 'G',
    'c': 'g',
}

COMPLEMENT_RNA = {
    'A': 'U',
    'a': 'u',
    'G': 'C',
    'g': 'c',

    'U': 'A',
    'u': 'a',
    'C': 'G',
    'c': 'g',
}


def reverse(sequences: list[str] | Tuple) -> list[str]:
    reverse_sequences = []
    for sequence in sequences:
        reverse_sequences.append(sequence[::-1])
    return reverse_sequences


def complement(sequences: list[str] | Tuple) -> list[str]:
    complemented = []
    for sequence in sequences:
        complemented.append('')
        if is_dna_or_rna(sequence):
            for nucl in sequence:
                complemented[-1] += COMPLEMENT_DNA[nucl]
        else:
            for nucl in sequence:
                complemented[-1] += COMPLEMENT_RNA[nucl]

    return complemented


def transcribe(sequences: list[str] | Tuple) -> list[str]:
    transcribe_sequences = []
    for sequence in sequences:
        transcribe_sequences.append('')
        for nucl in sequence:
            if nucl.upper() == 'A':
                transcribe_sequences[-1] += nucl
            elif nucl.upper() == 'G':
                transcribe_sequences[-1] += nucl
            elif nucl.upper() == 'C':
                transcribe_sequences[-1] += nucl
            elif nucl.upper() == 'T':
                if nucl == 't':
                    transcribe_sequences[-1] += 'u'
                elif nucl == 'T':
                    transcribe_sequences[-1] += 'U'

    return transcribe_sequences
