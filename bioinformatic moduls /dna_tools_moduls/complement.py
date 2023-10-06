def is_dna_rna(seq):
    for i in range(0, len(seq)):
        if seq[i].lower() == 't':
            return True
        elif seq[i].lower() == 'u':
            return False


complement_dna = {
    'A': 'T',
    'a': 't',
    'G': 'C',
    'g': 'c',

    'T': 'A',
    't': 'a',
    'C': 'G',
    'c': 'g',
}

complement_rna = {
    'A': 'U',
    'a': 'u',
    'G': 'C',
    'g': 'c',

    'U': 'A',
    'u': 'a',
    'C': 'G',
    'c': 'g',
}


def complement(sequences):
    complemented = []
    for i in range(0, len(sequences)):
        complemented.append('')
        if is_dna_rna(sequences[i]):
            for j in range(0, len(sequences[i])):
                complemented[i] += complement_dna[sequences[i][j]]
        else:
            for j in range(0, len(sequences[i])):
                complemented[i] += complement_rna[sequences[i][j]]

    return complemented
