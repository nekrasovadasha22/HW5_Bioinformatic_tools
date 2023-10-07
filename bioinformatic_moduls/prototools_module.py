"""
This is a prototool. WE ARE SORRY!!!
"""

from typing import List, Optional, Tuple, Union

AMINOACIDS_DICT = {
    'Ala': {'TO_1': 'A',
            'PROTEIN_TO_RNA_COMBINATION': {'GCU', 'GCC', 'GCA', 'GCG'},
            'PKA_AMINOACIDS': [2.34, 9.69],
            'MOLECULAR_WEIGHTS': 89},
    'Arg': {'TO_1': 'R',
            'PROTEIN_TO_RNA_COMBINATION': {'CGU', 'CGC', 'CGA', 'CGG', 'AGA',
                                           'AGG'},
            'PKA_AMINOACIDS': [2.17, 9.04, 12.68],
            'MOLECULAR_WEIGHTS': 174},
    'Asn': {'TO_1': 'N',
            'PROTEIN_TO_RNA_COMBINATION': {'AAU', 'AAC'},
            'PKA_AMINOACIDS': [1.88, 9.60, 3.65],
            'MOLECULAR_WEIGHTS': 132},
    'Asp': {'TO_1': 'D',
            'PROTEIN_TO_RNA_COMBINATION': {'GAU', 'GAC'},
            'PKA_AMINOACIDS': [1.88, 9.60, 3.65],
            'MOLECULAR_WEIGHTS': 133},
    'Cys': {'TO_1': 'C',
            'PROTEIN_TO_RNA_COMBINATION': {'UGU', 'UGC'},
            'PKA_AMINOACIDS': [1.96, 10.28, 8.18],
            'MOLECULAR_WEIGHTS': 121},
    'Glu': {'TO_1': 'Q',
            'PROTEIN_TO_RNA_COMBINATION': {'GAA', 'GAG'},
            'PKA_AMINOACIDS': [2.19, 9.67, 4.25],
            'MOLECULAR_WEIGHTS': 147},
    'Gln': {'TO_1': 'E',
            'PROTEIN_TO_RNA_COMBINATION': {'CAA', 'CAG'},
            'PKA_AMINOACIDS': [2.17, 9.13],
            'MOLECULAR_WEIGHTS': 146},
    'Gly': {'TO_1': 'G',
            'PROTEIN_TO_RNA_COMBINATION': {'GGU', 'GGC', 'GGA', 'GGG'},
            'PKA_AMINOACIDS': [2.34, 9.60],
            'MOLECULAR_WEIGHTS': 75},
    'His': {'TO_1': 'E',
            'PROTEIN_TO_RNA_COMBINATION': {'CAU', 'CAC'},
            'PKA_AMINOACIDS': [1.82, 9.17],
            'MOLECULAR_WEIGHTS': 155},
    'Ile': {'TO_1': 'I',
            'PROTEIN_TO_RNA_COMBINATION': {'AUU', 'AUC', 'AUA'},
            'PKA_AMINOACIDS': [2.36, 9.68],
            'MOLECULAR_WEIGHTS': 131},
    'Leu': {'TO_1': 'L',
            'PROTEIN_TO_RNA_COMBINATION': {'CUU', 'CUC', 'CUA', 'CUG'},
            'PKA_AMINOACIDS': [2.36, 9.60],
            'MOLECULAR_WEIGHTS': 131},
    'Lys': {'TO_1': 'K',
            'PROTEIN_TO_RNA_COMBINATION': {'AAA', 'AAG'},
            'PKA_AMINOACIDS': [2.18, 8.95, 10.53],
            'MOLECULAR_WEIGHTS': 146},
    'Met': {'TO_1': 'M',
            'PROTEIN_TO_RNA_COMBINATION': {'AUG'},
            'PKA_AMINOACIDS': [2.28, 9.21],
            'MOLECULAR_WEIGHTS': 149},
    'Phe': {'TO_1': 'F',
            'PROTEIN_TO_RNA_COMBINATION': {'UUU', 'UUC'},
            'PKA_AMINOACIDS': [2.20, 9.13],
            'MOLECULAR_WEIGHTS': 165},
    'Pro': {'TO_1': 'P',
            'PROTEIN_TO_RNA_COMBINATION': {'CCU', 'CCC', 'CCA', 'CCG'},
            'PKA_AMINOACIDS': [1.99, 10.96],
            'MOLECULAR_WEIGHTS': 115},
    'Ser': {'TO_1': 'S',
            'PROTEIN_TO_RNA_COMBINATION': {'UCU', 'UCC', 'UCA', 'UCG'},
            'PKA_AMINOACIDS': [2.21, 9.15],
            'MOLECULAR_WEIGHTS': 105},
    'Thr': {'TO_1': 'T',
            'PROTEIN_TO_RNA_COMBINATION': {'ACU', 'ACC', 'ACA', 'ACG'},
            'PKA_AMINOACIDS': [2.11, 9.62],
            'MOLECULAR_WEIGHTS': 119},
    'Tyr': {'TO_1': 'W',
            'PROTEIN_TO_RNA_COMBINATION': {'UAU', 'UAC'},
            'PKA_AMINOACIDS': [2.20, 9.11, 10.07],
            'MOLECULAR_WEIGHTS': 181},
    'Trp': {'TO_1': 'Y',
            'PROTEIN_TO_RNA_COMBINATION': {'UGG'},
            'PKA_AMINOACIDS': [2.38, 9.39],
            'MOLECULAR_WEIGHTS': 204},
    'Val': {'TO_1': 'V',
            'PROTEIN_TO_RNA_COMBINATION': {'GUU', 'GUC', 'GUA', 'GUG'},
            'PKA_AMINOACIDS': [2.32, 9.62],
            'MOLECULAR_WEIGHTS': 117},
}

# A dictionary where keys are 1-letter and values are 3-letters codes
TO_3_DICT = {nested_dict['TO_1']: key for key,
nested_dict in AMINOACIDS_DICT.items()}

TRANSCRIBE_DICT: dict = {'A': 'A',
                         'U': 'T',
                         'G': 'G',
                         'C': 'C',
                         'a': 'a',
                         'u': 't',
                         'g': 'g',
                         'c': 'c'}


def is_one_letter(seq: str) -> bool:
    """
    Defines whether the sequence is 1 coded.

    Args:
    - seq - sequence to check

    Returns:
    - bool
    """
    return all(aa.isalpha() and aa.isupper() for aa in seq)


def recode(seq: str) -> str:
    """
    Translate 1-letter to 3-letter encoding if 1-letter
    encoded sequence is given and vice versa.

    Args:
    - seq - sequence or list of sequences to recode

    Returns:
    - function_result - a dictionary containing recoded sequences as values
    for original sequences keys
    """

    if is_one_letter(seq):
        # Translate 1-letter to 3-letter coded sequence
        three_letter_sequence = ""
        for aa in seq:
            three_letter_code = TO_3_DICT.get(aa, aa)
            three_letter_sequence += three_letter_code
        return three_letter_sequence
    # Translate 3-letter to 1-letter coded sequence
    one_letter_sequence = ""
    for aa in range(0, len(seq), 3):
        amino_acid = seq[aa:aa + 3]
        one_letter_sequence += AMINOACIDS_DICT[amino_acid]['TO_1']
    return one_letter_sequence


def prettify_alignment(aligned_seq_on: str, aligned_seq2: str) -> None:
    """
    Prettifies alignment output by printing out two
    sequences on top of each other

    Finds the start of aligned sequence in the longer of sequences.\\
    Prints the longer sequence as an upper one and aligned sequence
    is bellow separated via vertical lines

    Args:
    - aligned_seq_on, aligned_seq2 - sequences
    from the local_alignment()

    Returns:
    None \\
    Prints out the prettified view in stdout
    """

    print(aligned_seq_on)
    print('|' * len(aligned_seq2))
    print(aligned_seq2)


def count_protein_molecular_weight(*seqs_list: Union[List[str], str]) -> dict:
    """
    :param seqs_list: seqs_list is a list of strings without whitespace
    (e.g. 'AlaSer'). You can put as many sequences as you wish.
    :return: This function returns molecular weight of the protein.
    """
    result = {}
    for seq in seqs_list:
        protein_weight = 0
        aminoacids = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        for i, aminoacid in enumerate(aminoacids):
            if aminoacid in AMINOACIDS_DICT.keys():
                aminoacid_weight = (AMINOACIDS_DICT[aminoacid]
                ['MOLECULAR_WEIGHTS'])
                protein_weight += aminoacid_weight
                result[seq] = protein_weight
    return result


def from_proteins_seqs_to_rna(*seqs_list: Union[List[str], str]) -> dict:
    """
    :param seqs_list: a list of strings with type 'ValTyrAla','AsnAspCys'.
    You can pass more than one sequence at the time.
    :return: dictionary, where [key] is your input protein sequences
    and values are combinations of RNA codones, which encode this protein
    """

    answer_dictionary = {}
    for seq in seqs_list:

        rna_combination = ''
        divided_acids = [seq[i:i + 3] for i in range(0,
                                                     len(seq),
                                                     3)]
        for divided_acid in divided_acids:

            if divided_acid in AMINOACIDS_DICT.keys():
                rna_combination += next(iter(AMINOACIDS_DICT[divided_acid]
                                             ['PROTEIN_TO_RNA_COMBINATION']))
            else:
                raise ValueError('Non-protein aminoacids in sequence')
        answer_dictionary[seq] = rna_combination
    return answer_dictionary


def isoelectric_point_determination(*seqs_list: List[str] | str) -> dict:
    """
    :param seqs_list: a list of strings with type 'ValTyrAla','AsnAspCys'.
    You can pass more than one sequence at a time.
    :return: dictionary, where [key] is your input protein sequence and value
    is an isoelectric point of your input proteins
    """
    answer_dictionary = {}

    for aminoacids in seqs_list:
        divided_acids = [aminoacids[i:i + 3] for i in range(0,
                                                            len(aminoacids),
                                                            3)]
        for divided_acid in divided_acids:
            if divided_acid not in AMINOACIDS_DICT.keys():
                raise ValueError('Non-protein aminoacids in sequence')

        isoelectric_point_mean = 0
        count_groups = 0
        for acid_index, aminoacid in enumerate(divided_acids):
            if acid_index == 0:
                isoelectric_point_mean \
                    += (AMINOACIDS_DICT[aminoacid]['PKA_AMINOACIDS'][0])
                count_groups += 1
            elif acid_index == len(divided_acids) - 1:
                isoelectric_point_mean = (isoelectric_point_mean
                                          + (AMINOACIDS_DICT[aminoacid]
                        ['PKA_AMINOACIDS'][-1]))
                count_groups += 1
            else:
                if len(AMINOACIDS_DICT[aminoacid]['PKA_AMINOACIDS']) > 2:
                    isoelectric_point_mean = (isoelectric_point_mean
                                              + (AMINOACIDS_DICT[aminoacid]
                            ['PKA_AMINOACIDS'][1]))
                    count_groups += 1
        answer_dictionary[aminoacids] = isoelectric_point_mean / count_groups
    return answer_dictionary


def back_transcribe(*seqs_list: Union[List[str], str]) -> dict:
    """
    :param seqs_list: is a list of strings without whitespace. 
    You can put as many sequences as you wish.
    :return: This function returns a dictonary where key is inputed protein
    sequence and values are DNA codons
    """
    result = {}
    for seq in seqs_list:
        rna = list((from_proteins_seqs_to_rna(seq)).get(seq))
        for i in range(len(rna)):
            if rna[i] in TRANSCRIBE_DICT.keys():
                rna[i] = TRANSCRIBE_DICT[rna[i]]
        result[seq] = "".join(rna)
    return result


def count_gc_content(*seqs_list: Union[List[str], str]) -> dict:
    """
    :param seqs_list: is a list of strings without whitespace. 
    You can put as many sequences as you wish.
    :return: This function returns GC-content of DNA sequence, which encodes
    the protein
    """
    result = {}
    for seq in seqs_list:
        dna = list((back_transcribe(seq)).get(seq))
        gc_content = round(100 * (dna.count('G') + dna.count('C')) / len(dna))
        result[seq] = gc_content
    return result


def check_input(*args: List[str] | str, method: str) -> \
        Tuple[List[str], Optional[str]]:
    """
    Function to check the validity of the input.

    Args:
    - *args - are supposed to be all sequences to process
    - method - the method to process with method

    Returns:
    - seqs_list - list of sequences
    - seq_on (optional) - in case of local_alignment method
    """

    if len(args) == 0:
        # Handle the case where there are no arguments
        raise ValueError('No input defined.')
    else:
        if method not in ['recode',
                          'local_alignment',
                          'from_proteins_seqs_to_rna',
                          'isoelectric_point_determination',
                          'count_protein_molecular_weight',
                          'back_transcribe',
                          'count_gc_content']:
            raise ValueError(method, ' is not a valid method.')
        else:
            # Form a list with sequences from the input
            seqs_list = list(args)
            if method == 'local_alignment':
                if len(seqs_list) < 2:
                    raise IndexError('Need at least two sequences to align.')
                for i, seq in enumerate(seqs_list):
                    if not is_one_letter(seq):
                        print('Warning! Function local_alignment() needs '
                              '1-letter encoded sequences. Your sequence '
                              'will be mutated to a 1-letter encoding.')
                        seqs_list[i] = recode(seq)
                        print(seq, ' sequence has been mutated into: ',
                              seqs_list[i])
                        seq_on = seqs_list.pop(0)
                return seqs_list, seq_on
            for i, seq in enumerate(seqs_list):
                if is_one_letter(seq):
                    print(f'Warning! Function {method}() needs '
                          '3-letter encoded sequences. Your sequence '
                          'will be mutated to a 3-letter encoding.')
                    seqs_list[i] = recode(seq)
                    print(seq, ' sequence has been mutated into: ',
                          seqs_list[i])
            seq_on = None
            return seqs_list, seq_on


def main(*args: Tuple[Union[List[str], str]],
         method: Optional[str] = None) -> dict:
    """
    This function provides the access to the following methods:

    1. Translate 1 letter to 3 letter encoding and vice versa - the last
    argument: 'recode'
        - needs at least 1 sequence 1- or 3- letter encoded. Can recive
        more than 1 sequences
        - returns a dictionary containing translations between 1- and 3-
        letter codes
    2. Find possible RNA sequences for defined protein sequence - the
    last argument: 'from_proteins_seqs_to_rna'
        - needs at least 1 protein sequence 3-letter encoded
        - returns a dictionary, where key is your input protein sequences
        and values are combinations of RNA codones, which encode this protein

    3. Determinate isoelectric point - the last argument:
    'isoelectric_point_determination'
        - needs an input containing at least 1 aminoacid. Can recive multiple
        different protein sequences
        - returns a dictionary, where key is your input protein sequence and
        value is an isoelectric point of this protein

    4. Calculate protein molecular weight - the last argument:
    'count_protein_molecular_weight'
        - Seqs is an argument of the function. It is a string without
    whitespace (e.g. 'AlaSer'). You can put as many arguments as you wish.
        - returns a dictionary with protein sequences as keys and their
        calculated molecular weight as corresponding values

    5. Determine possible DNA sequence from protein sequence - the last
    argument: 'back_transcribe'
        - needs a string without whitespaces. You can put as many arguments as
        you wish.
        - returns a dictonary where keys are inputed protein sequences and
        corresponding values are possible DNA codons

    6. Calculate a GC ratio in a possible DNA sequence of a given aminoacid
    sequence - the last argument 'count_gc_content'
        - needs a string without whitespaces. You can put as many sequences
        as you wish.
        - returns a dictionary where keys are inputed aminoacid sequences and
        GC-content of DNA sequence, which encodes the protein are
        corresponding values

    Args:
    - *args - are supposed to be all sequences to process
    - method is a kwarg - the method to process with.

    Returns:
    function_result - a dictionary with the result of a chosen function
    """

    seqs_list, seq_on = check_input(*args, method=method)
    print(f'Your sequences are: {seqs_list}',
          f'The method is: {method}', sep='\n')

    match method:

        case 'recode':

            recode_dict: dict = {}
            for seq in seqs_list:
                recode_dict[seq] = recode(seq=seq)
            return recode_dict

        case 'from_proteins_seqs_to_rna':

            return from_proteins_seqs_to_rna(*seqs_list)

        case 'count_protein_molecular_weight':

            return count_protein_molecular_weight(*seqs_list)

        case 'isoelectric_point_determination':

            return isoelectric_point_determination(*seqs_list)

        case 'back_transcribe':

            return back_transcribe(*seqs_list)

        case 'count_gc_content':

            return count_gc_content(*seqs_list)