def transcribe(sequences):
    transcribe_sequences = []
    for i in range(0, len(sequences)):
        transcribe_sequences.append('')
        for j in range(0, len(sequences[i])):
            if sequences[i][j].upper() == 'A':
                transcribe_sequences[i] += sequences[i][j]
            elif sequences[i][j].upper() == 'G':
                transcribe_sequences[i] += sequences[i][j]
            elif sequences[i][j].upper() == 'C':
                transcribe_sequences[i] += sequences[i][j]
            elif sequences[i][j].upper() == 'T':
                if sequences[i][j] == 't':
                    transcribe_sequences[i] += 'u'
                elif sequences[i][j] == 'T':
                    transcribe_sequences[i] += 'U'

    return transcribe_sequences
