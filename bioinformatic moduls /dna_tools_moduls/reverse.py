def reverse(sequences):
    reverse_sequences = []
    for i in range(0, len(sequences)):
        reverse_sequences.append('')
        for j in range(len(sequences[i]) - 1, -1, -1):
            reverse_sequences[i] += sequences[i][j]
    return reverse_sequences
