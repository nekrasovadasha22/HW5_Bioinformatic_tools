import os

OUTPUT_FASTA_SUFFIX = '_output_file'
FASTA_FILE_EXT = '.fasta'

LOCUS_GBK_KEYWORD = 'LOCUS'
FEATURES_GBK_KEYWORD = 'FEATURES'
ORIGIN_GBK_KEYWORD = 'ORIGIN'
CDS_GBK_KEYWORD = 'CDS'
END_OF_PART_GBK_KEYWORD = '//'

GENE_CDS_KEYWORD = '/gene'
TRANSLATION_CDS_KEYWORD = '/translation'

def get_last_dict_value(input_dict: dict):
    keys = list(input_dict.keys())
    return input_dict[keys[-1]]


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = ''):
    """
       Converts a multi-line FASTA file to a one-line FASTA format.

       Args:
           input_fasta (str): The path to the input multi-line FASTA file.
           output_fasta (str, optional): The path to the output one-line FASTA
           file. If not provided, it will be generated based on the input file.

       Raises:
           ValueError: If the input FASTA file does not exist.

       Note:
           The function reads the input multi-line FASTA file and converts it
           to a one-line FASTA file, where each sequence is on a single line.

       """
    if not os.path.exists(input_fasta):
        raise ValueError(f'Your {input_fasta} doesn\'t exist')
    if output_fasta != '':
        if not output_fasta.endswith(FASTA_FILE_EXT):
            output_fasta += FASTA_FILE_EXT
    else:
        output_fasta = os.path.join(os.path.dirname(input_fasta),
                                    os.path.splitext
                                    (input_fasta)[0] +
                                    OUTPUT_FASTA_SUFFIX + FASTA_FILE_EXT)
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


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list[str],
                                   n_before: int = 1,
                                   n_after: int = 1, output_fasta: str = ''):
    """
       Selects genes from a GenBank file and writes them to a FASTA file.

       Args:
           input_gbk (str): The path to the input GenBank file.
           genes (list[str]): A list of gene names to be selected.
           n_before (int, optional): The number of genes before the selected
           gene.
           n_after (int, optional): The number of genes after the selected
           gene.
           output_fasta (str, optional): The path to the output FASTA file.
               If not provided, it will be generated based on the input file.

       Raises:
           ValueError: If the input GenBank file does not exist, genes list is
           empty,
           or if n_before or n_after is non-positive.

       Note:
           The function reads the input GenBank file, selects the specified
           genes along with a certain number of genes before and after them,
           and writes the selected genes to a FASTA file.
       """
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
                                    OUTPUT_FASTA_SUFFIX + FASTA_FILE_EXT)
    parsed_data = {}
    with (open(input_gbk) as input_gbk_file):
        current_line = input_gbk_file.readline()
        while current_line != '':

            if current_line.startswith(LOCUS_GBK_KEYWORD):
                parsed_data[current_line.split()[1]] = ['', '', []]

            while not current_line.startswith(FEATURES_GBK_KEYWORD):
                current_line = input_gbk_file.readline()

            current_line = input_gbk_file.readline()
            separate_line = current_line.split()
            get_last_dict_value(parsed_data)[0] = separate_line[1]
            while not current_line.startswith(ORIGIN_GBK_KEYWORD):
                current_line = input_gbk_file.readline()
                separate_line = current_line.split()

                if separate_line[0] == CDS_GBK_KEYWORD:
                    get_last_dict_value(parsed_data)[2] \
                        .append([separate_line[1], '', ''])

                    current_line = input_gbk_file.readline().strip()
                    separate_line = current_line.split()
                    while separate_line[0] != CDS_GBK_KEYWORD:
                        cds_parameter = separate_line[0].split('=')
                        if cds_parameter[0] == GENE_CDS_KEYWORD:
                            get_last_dict_value(parsed_data)[2][-1][1] \
                                = cds_parameter[1].replace('"', '')
                        if cds_parameter[0] == TRANSLATION_CDS_KEYWORD:
                            translation = [cds_parameter[1]]
                            while not current_line.endswith('"'):
                                current_line = input_gbk_file.readline().strip()
                                separate_line = current_line.split()
                                translation.append(separate_line[0])
                            translation_line = ''.join(translation)
                            get_last_dict_value(parsed_data)[2][-1][2] \
                                = translation_line.replace('"', '')
                            break
                        current_line = input_gbk_file.readline()
                        separate_line = current_line.split()

            current_line = input_gbk_file.readline().strip()
            while current_line != END_OF_PART_GBK_KEYWORD:
                current_line = input_gbk_file.readline().strip()
            current_line = input_gbk_file.readline()

            with open(output_fasta, mode="w") as output_fasta_file:
                for current_key in parsed_data:
                    for i in range(0, len(parsed_data[current_key][2])):
                        if any(parsed_data[current_key][2][i][1] ==
                                gene for gene in genes):
                            left_translation = []
                            right_translation = []
                            for k in range(i - 1, i - n_before - 1, -1):
                                if k < 0:
                                    break
                                left_translation\
                                    .append(parsed_data[current_key][2][k][2])
                            for j in range(i + 1, i + n_after + 1):
                                if i >= len(parsed_data[current_key][2]):
                                    break
                                right_translation\
                                    .append(parsed_data[current_key][2][j][2])
                            left_join = ''.join(left_translation[::-1])
                            right_join = ''.join(right_translation)
                            merged_translation = left_join + right_join
                            output_fasta_file.write(current_key + '\n'\
                                                   + merged_translation + '\n')


