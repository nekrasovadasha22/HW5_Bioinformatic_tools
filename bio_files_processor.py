import os
from typing import Tuple


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


def select_genes_from_gbk_to_fasta(input_gbk:str, genes:str, n_before:int = 1,
                                   n_after:int = 1, output_fasta:str=''):
    if not os.path.exists(input_gbk):
        raise ValueError(f'Your {input_gbk} doesn\'t exist')
    if genes == '':
        raise ValueError(f'{genes} is empty!')
    if n_before <= 0 or n_after <= 0:
        raise ValueError('Enter the genes quantity before and after \
        interesting gene')
    if output_fasta == '':
        output_fasta = os.path.join(os.path.dirname(input_gbk),
                                    os.path.splitext
                                    (input_gbk)[0] +
                                    '_output_file' + '.fasta')
    parsed_data = {}
    with open(input_gbk) as input_gbk_file:
        current_line = input_gbk_file.readline()
        while current_line != '':
            if current_line.startswith('LOCUS'):
                parsed_data[current_line.split()[1]] = ('', '', [])
            while not current_line.startswith('FEATURES'):
                current_line = input_gbk_file.readline()
            current_line = input_gbk_file.readline()
            separate_line = current_line.split()
            parsed_data[parsed_data.keys()[-1]][0] = separate_line[1]
            while True and current_line != '//':
                current_line = input_gbk_file.readline()
                separate_line = current_line.split()





