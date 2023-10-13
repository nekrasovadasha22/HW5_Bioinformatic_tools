import os


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = ''):
    if not os.path.exists(input_fasta):
        raise ValueError(f'Your {input_fasta} doesn\'t exist')
    if output_fasta != '':
        if not output_fasta.endswith('.fasta'):
            output_fasta += '.fasta'
    else:
        output_fasta = os.path.join(os.path.dirname(input_fasta),
                                    os.path.splitext
                                    (input_fasta)[0] +
                                    '_output_file' + '.fasta')
    with open(input_fasta) as input_fasta_file:
        with open(output_fasta, mode="w") as output_fasta_file:
            current_line = input_fasta_file.readline().strip()
            while current_line != '':
                output_fasta_file.write(f'{current_line}\n')
                sequence_lines = []
                current_line = input_fasta_file.readline().strip()
                sequence_lines.append(current_line)
                while not current_line.startswith('>') and current_line != '':
                    current_line = input_fasta_file.readline().strip()
                    if not current_line.startswith('>') and current_line != '':
                        sequence_lines.append(current_line)
                sequence_lines = ''.join(sequence_lines)
                output_fasta_file.write(f'{sequence_lines}\n')


