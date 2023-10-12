import os.path
from typing import Tuple, Optional, Union, List
from bioinformatic_moduls.dna_tools_module import (reverse,
                                                   transcribe, complement)
from bioinformatic_moduls.fastqtools_module import (phread33_converter,
                                                    convert_fastq_to_dict)
from bioinformatic_moduls.prototools_module \
    import (convert_aa_coding, get_protein_mrnas,
            calc_protein_molecular_weight,
            isoelectric_point_calculating, back_transcribe, count_gc_content,
            check_input)


def run_dna_rna_tools(*parameters: str) -> list[str] | str:
    """
        Run DNA and RNA sequence manipulation tools.

        This function accepts a variable number of parameters, with the last parameter
        specifying the name of the tool to be used. The preceding parameters should
        contain one or more DNA or RNA sequences as strings, represented by a combination
        of characters from the set ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u'].

        Parameters:
        *parameters (str): Variable number of DNA or RNA sequences and the tool name.

        Returns:
        [list[str], str]: Depending on the tool used, it returns a list of
        sequences or a single sequence as a string. If only one sequence is processed,
        it returns the sequence as a string. If the tool name is invalid or the answer
        is None, it raises a ValueError.

        Raises:
        - ValueError: If the parameters are None, there are not enough parameters,
          a given parameter is None, the sequences contain characters other than
          ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u'], or the tool name is unknown.
        - ValueError: If an RNA sequence contains 'T' or a DNA sequence contains 'U'.
        - ValueError: If the answer is None or does not contain valid sequences.

        Example:
        >>> run_dna_rna_tools("ATGC", "AUG", "transcribe")
        ['AUGC', 'AUG', 'transcribe']

        >>> run_dna_rna_tools("ATGC", "UAGC", "transcribe")
        Traceback (most recent call last):
          ...
        ValueError: RNA sequence cannot contain T
        """
    available_rna_dna_symbols = \
        ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c', 'U', 'u']

    if parameters is None:
        raise ValueError('Parameters are None!')
    if len(parameters) < 2:
        raise ValueError('Parameters are not enough!')

    tool_name = parameters[-1]
    sequences = parameters[:-1]

    for sequence in sequences:
        if sequence is None:
            raise ValueError('Given parameter was None')
        else:
            for nucl in sequence:
                if not (available_rna_dna_symbols
                        .__contains__(nucl)):
                    raise ValueError('Parameters are not nucleotide sequences')

    for sequence in sequences:
        is_dna = False
        is_rna = False
        for nucl in sequence:
            if nucl.upper() == 'T':
                if not is_rna:
                    is_dna = True
                else:
                    raise ValueError('RNA sequence cannot contain T')
            if nucl.upper() == 'U':
                if not is_dna:
                    is_rna = True
                else:
                    raise ValueError('DNA sequence cannot contain U')
    answer = None

    if tool_name == 'transcribe':
        answer = transcribe(sequences)
    elif tool_name == 'reverse':
        answer = reverse(sequences)
    elif tool_name == 'complement':
        answer = complement(sequences)
    elif tool_name == 'reverse_complement':
        complemented = complement(sequences)
        answer = reverse(complemented)
    else:
        raise ValueError('Unknown tool')

    if len(answer) < 2:
        return answer[0]
    elif answer is None:
        raise ValueError('Answer is None')
    else:
        return answer


def run_fastq_filter(input_path:str, gc_bounds: int | float | Tuple = (20, 80),
                     length_bounds: int | float | Tuple = (0, 2 ** 32),
                     quality_threshold: int = 0, output_filename = '') -> dict:
    """
 Filter sequences from a FastQ file based on various criteria.

 This function filters sequences from a FastQ file represented as a dictionary
 with sequence names as keys and values as lists containing the sequence data,
 quality scores. Sequences are filtered based on GC content, sequence length,
 and quality score.

 Parameters:
 seqs (dict): A dictionary containing FastQ sequences, where keys are sequence
     names, and values are lists with at least two elements, where the first
     element is the sequence string and the last element is the quality score
     string.
 gc_bounds (int | float | Tuple, optional): The GC content bounds for filtering.
     Default is (20, 80), indicating sequences with GC content between 20% and 80%
     (inclusive) are retained. Can be a single value (lower bound) or a tuple
     (lower and upper bounds).
 length_bounds (int | float | Tuple, optional): The sequence length bounds for
     filtering. Default is (0, 2^32), indicating no length filtering. Can be a
     single value (lower bound) or a tuple (lower and upper bounds).
 quality_threshold (int): The minimum quality score threshold for
     filtering sequences. Default is 0, meaning all sequences pass this filter.

 Returns:
 dict: A dictionary containing filtered FastQ sequences, preserving the original
 keys.

 Raises:
 - ValueError: If seqs is None or empty.
 - ValueError: If the quality score string for a sequence has a length of 0.
 """
    seqs = convert_fastq_to_dict(input_path)
    if seqs is None:
        raise ValueError('Your fastq_files are None')
    elif len(seqs) == 0:
        raise ValueError('Your fastq_files are empty')

    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)

    fastq_filtered = {}
    for fastq_name in seqs:
        encrypted_phread_score = seqs[fastq_name][-1]
        if len(encrypted_phread_score) == 0:
            raise ValueError(f'In your fastq {fastq_name} '
                             f'file pherad score string is length 0')
        sequence = seqs[fastq_name][0]

        converted_phread_score = phread33_converter(encrypted_phread_score)
        if converted_phread_score >= quality_threshold and \
                length_bounds[0] <= len(sequence) <= \
                length_bounds[1]:
            gc_content = ((sequence.count('G')
                           + sequence.count('C'))
                          / len(sequence) * 100)
            if gc_bounds[1] >= gc_content >= gc_bounds[0]:
                fastq_filtered[fastq_name] = seqs[fastq_name]

    if output_filename == '':
        output_filename = os.path.basename(input_path)

    if '.fastq' in output_filename:
        output_filename = output_filename.replace('.fastq', '')

    output_dir = os.path.join(os.path.dirname(input_path), 'fastq_filtrator_results')

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(os.path.join(output_dir, f'{output_filename}.fastq'), mode='w') as result_file:
        for seq_name in fastq_filtered:
            result_file.write(f'{seq_name}\n')
            result_file.write(f'{fastq_filtered[seq_name][0]}\n')
            q_score_name = seq_name.replace('@', '+')
            result_file.write(f'{q_score_name}\n')
            result_file.write(f'{fastq_filtered[seq_name][1]}\n')



def run_prototools(*args: List[str] | str,
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

        case 'convert_aa_coding':

            recode_dict: dict = {}
            for seq in seqs_list:
                recode_dict[seq] = convert_aa_coding(seq=seq)
            return recode_dict

        case 'get_protein_mrnas':

            return get_protein_mrnas(*seqs_list)

        case 'calc_protein_molecular_weight':

            return calc_protein_molecular_weight(*seqs_list)

        case 'isoelectric_point_calculating':

            return isoelectric_point_calculating(*seqs_list)

        case 'back_transcribe':

            return back_transcribe(*seqs_list)

        case 'count_gc_content':

            return count_gc_content(*seqs_list)


resutl = run_fastq_filter("/home/daria/repos/.fr-XATXC2/HW6_Files-main/example_data/example_fastq.fastq")
print(resutl)