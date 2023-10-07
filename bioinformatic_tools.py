from typing import Tuple, Optional, Union, List
from bioinformatic_moduls.dna_tools_module import (reverse,
                                                   transcribe, complement)
from bioinformatic_moduls.fastqtools_module import phread33_converter
from bioinformatic_moduls.prototools_module \
    import (recode, from_proteins_seqs_to_rna, count_protein_molecular_weight,
            isoelectric_point_determination, back_transcribe, count_gc_content,
            check_input)


def run_dna_rna_tools(*parameters: str) -> list[str] | str:
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


def run_fastq_filter(seqs: dict, gc_bounds: int | float | Tuple = (20, 80),
                     length_bounds: int | float | Tuple = (0, 2 ** 32),
                     quality_threshold: int = 0) -> dict:
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
    return fastq_filtered


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
