class NucleicAcidTools:
    def is_dna_or_rna(self, seq: str) -> bool:
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

    def reverse(self, sequences: list[str] | Tuple) -> list[str]:
        reverse_sequences = []
        for sequence in sequences:
            reverse_sequences.append(sequence[::-1])
        return reverse_sequences

    def complement(self, sequences: list[str] | Tuple) -> list[str]:
        complemented = []
        for sequence in sequences:
            complemented.append('')
            if self.is_dna_or_rna(sequence):
                for nucl in sequence:
                    complemented[-1] += self.COMPLEMENT_DNA[nucl]
            else:
                for nucl in sequence:
                    complemented[-1] += self.COMPLEMENT_RNA[nucl]

        return complemented

    def transcribe(self, sequences: list[str] | Tuple) -> list[str]:
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


class FastqTools:
    def phread33_converter(self, seqs: str) -> float | int:
        phread_score = 0
        for nucl_quality in seqs:
            phread_score += ord(nucl_quality) - 33
        return phread_score / len(seqs)

    def convert_fastq_to_dict(self, input_path: str) -> dict:
        seqs = {}
        fastq_file_lines = []
        with open(input_path) as fastq_file:
            fastq_file_lines = fastq_file.readlines()
        for i in range(0, len(fastq_file_lines), 4):
            seqs[fastq_file_lines[i].strip()] = (
                fastq_file_lines[i + 1].strip(),
                fastq_file_lines[i + 3].strip())
        return seqs


class ProteinTools:
    AMINOACIDS_DICT = {
        'Ala': {'TO_1': 'A',
                'PROTEIN_TO_RNA_COMBINATION': {'GCU', 'GCC', 'GCA', 'GCG'},
                'PKA_AMINOACIDS': [2.34, 9.69],
                'MOLECULAR_WEIGHTS': 89},
        'Arg': {'TO_1': 'R',
                'PROTEIN_TO_RNA_COMBINATION': {'CGU', 'CGC', 'CGA', 'CGG',
                                               'AGA',
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

    def is_one_letter(self, seq: str) -> bool:
        """
        Defines whether the sequence is 1 coded.

        Args:
        - seq - sequence to check

        Returns:
        - bool
        """
        return all(aa.isalpha() and aa.isupper() for aa in seq)

    def convert_aa_coding(self, seq: str) -> str:
        """
        Translate 1-letter to 3-letter encoding if 1-letter
        encoded sequence is given and vice versa.

        Args:
        - seq - sequence or list of sequences to recode

        Returns:
        - function_result - a dictionary containing recoded sequences as values
        for original sequences keys
        """

        if self.is_one_letter(seq):
            three_letter_sequence = ""
            for aa in seq:
                three_letter_code = self.TO_3_DICT.get(aa, aa)
                three_letter_sequence += three_letter_code
            return three_letter_sequence
        one_letter_sequence = ""
        for aa in range(0, len(seq), 3):
            amino_acid = seq[aa:aa + 3]
            one_letter_sequence += self.AMINOACIDS_DICT[amino_acid]['TO_1']
        return one_letter_sequence

    def calc_protein_molecular_weight(self,
                                      *seqs_list: Union[
                                          List[str], str]) -> dict:
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
                if aminoacid in self.AMINOACIDS_DICT.keys():
                    aminoacid_weight = (self.AMINOACIDS_DICT[aminoacid]
                    ['MOLECULAR_WEIGHTS'])
                    protein_weight += aminoacid_weight
                    result[seq] = protein_weight
        return result

    def get_protein_mrnas(self, *seqs_list: Union[List[str], str]) -> dict:
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

                if divided_acid in self.AMINOACIDS_DICT.keys():
                    rna_combination += next(
                        iter(self.AMINOACIDS_DICT[divided_acid]
                             [
                                 'PROTEIN_TO_RNA_COMBINATION']))
                else:
                    raise ValueError('Non-protein aminoacids in sequence')
            answer_dictionary[seq] = rna_combination
        return answer_dictionary

    def isoelectric_point_calculating(self,
                                      *seqs_list: List[str] | str) -> dict:
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
                if divided_acid not in self.AMINOACIDS_DICT.keys():
                    raise ValueError('Non-protein aminoacids in sequence')

            isoelectric_point = 0
            count_groups = 0
            for acid_index, aminoacid in enumerate(divided_acids):
                if acid_index == 0:
                    isoelectric_point \
                        += (
                    self.AMINOACIDS_DICT[aminoacid]['PKA_AMINOACIDS'][0])
                    count_groups += 1
                elif acid_index == len(divided_acids) - 1:
                    isoelectric_point = (isoelectric_point
                                         + (self.AMINOACIDS_DICT[aminoacid]
                            ['PKA_AMINOACIDS'][-1]))
                    count_groups += 1
                else:
                    if len(self.AMINOACIDS_DICT[aminoacid][
                               'PKA_AMINOACIDS']) > 2:
                        isoelectric_point = (isoelectric_point
                                             + (self.AMINOACIDS_DICT[aminoacid]
                                ['PKA_AMINOACIDS'][1]))
                        count_groups += 1
            answer_dictionary[aminoacids] = isoelectric_point / count_groups
        return answer_dictionary

    def back_transcribe(self, *seqs_list: Union[List[str], str]) -> dict:
        """
        :param seqs_list: is a list of strings without whitespace.
        You can put as many sequences as you wish.
        :return: This function returns a dictonary where key is inputed protein
        sequence and values are DNA codons
        """
        result = {}
        for seq in seqs_list:
            rna = list((self.get_protein_mrnas(seq)).get(seq))
            for i in range(len(rna)):
                if rna[i] in self.TRANSCRIBE_DICT.keys():
                    rna[i] = self.TRANSCRIBE_DICT[rna[i]]
            result[seq] = "".join(rna)
        return result

    def count_gc_content(self, *seqs_list: Union[List[str], str]) -> dict:
        """
        :param seqs_list: is a list of strings without whitespace.
        You can put as many sequences as you wish.
        :return: This function returns GC-content of DNA sequence, which encodes
        the protein
        """
        result = {}
        for seq in seqs_list:
            dna = list((self.back_transcribe(seq)).get(seq))
            gc_content = round(
                100 * (dna.count('G') + dna.count('C')) / len(dna))
            result[seq] = gc_content
        return result

    def check_input(self, *args: Union[List[str], str], method: str) -> \
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
            if method not in ['convert_aa_coding',
                              'get_protein_mrnas',
                              'isoelectric_point_calculating',
                              'calc_protein_molecular_weight',
                              'back_transcribe',
                              'count_gc_content']:
                raise ValueError(method, ' is not a valid method.')
            else:
                # Form a list with sequences from the input
                seqs_list = list(args)
                for i, seq in enumerate(seqs_list):
                    if self.is_one_letter(seq):
                        print(f'Warning! Function {method}() needs '
                              '3-letter encoded sequences. Your sequence '
                              'will be mutated to a 3-letter encoding.')
                        seqs_list[i] = self.convert_aa_coding(seq)
                        print(seq, ' sequence has been mutated into: ',
                              seqs_list[i])
                        seq_on = None
                return seqs_list, seq_on


class FastqParser:
    OUTPUT_FASTA_SUFFIX = '_output_file'
    FASTA_FILE_EXT = '.fasta'

    LOCUS_GBK_KEYWORD = 'LOCUS'
    FEATURES_GBK_KEYWORD = 'FEATURES'
    ORIGIN_GBK_KEYWORD = 'ORIGIN'
    CDS_GBK_KEYWORD = 'CDS'
    END_OF_PART_GBK_KEYWORD = '//'

    GENE_CDS_KEYWORD = '/gene'
    TRANSLATION_CDS_KEYWORD = '/translation'

    def get_last_dict_value(self, input_dict: dict):
        keys = list(input_dict.keys())
        return input_dict[keys[-1]]

    def convert_multiline_fasta_to_oneline(self, input_fasta: str,
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
            if not output_fasta.endswith(self.FASTA_FILE_EXT):
                output_fasta += self.FASTA_FILE_EXT
        else:
            output_fasta = os.path.join(os.path.dirname(input_fasta),
                                        os.path.splitext
                                        (input_fasta)[0] +
                                        self.OUTPUT_FASTA_SUFFIX + self.FASTA_FILE_EXT)
        with open(input_fasta) as input_fasta_file:
            with open(output_fasta, mode="w") as output_fasta_file:
                current_line = input_fasta_file.readline().strip()
                while current_line != '':
                    output_fasta_file.write(f'{current_line}\n')
                    sequence_lines = []
                    current_line = input_fasta_file.readline().strip()
                    sequence_lines.append(current_line)
                    while not current_line.startswith(
                            '>') and current_line != '':
                        current_line = input_fasta_file.readline().strip()
                        if not current_line.startswith(
                                '>') and current_line != '':
                            sequence_lines.append(current_line)
                    sequence_lines = ''.join(sequence_lines)
                    output_fasta_file.write(f'{sequence_lines}\n')

    def select_genes_from_gbk_to_fasta(self, input_gbk: str, genes: list[str],
                                       n_before: int = 1,
                                       n_after: int = 1,
                                       output_fasta: str = ''):
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
                                        self.OUTPUT_FASTA_SUFFIX + self.FASTA_FILE_EXT)
        parsed_data = {}
        with (open(input_gbk) as input_gbk_file):
            current_line = input_gbk_file.readline()
            while current_line != '':

                if current_line.startswith(self.LOCUS_GBK_KEYWORD):
                    parsed_data[current_line.split()[1]] = ['', '', []]

                while not current_line.startswith(self.FEATURES_GBK_KEYWORD):
                    current_line = input_gbk_file.readline()

                current_line = input_gbk_file.readline()
                separate_line = current_line.split()
                self.get_last_dict_value(parsed_data)[0] = separate_line[1]
                while not current_line.startswith(self.ORIGIN_GBK_KEYWORD):
                    current_line = input_gbk_file.readline()
                    separate_line = current_line.split()

                    if separate_line[0] == self.CDS_GBK_KEYWORD:
                        self.get_last_dict_value(parsed_data)[2] \
                            .append([separate_line[1], '', ''])

                        current_line = input_gbk_file.readline().strip()
                        separate_line = current_line.split()
                        while separate_line[0] != self.CDS_GBK_KEYWORD:
                            cds_parameter = separate_line[0].split('=')
                            if cds_parameter[0] == self.GENE_CDS_KEYWORD:
                                self.get_last_dict_value(parsed_data)[2][-1][1] \
                                    = cds_parameter[1].replace('"', '')
                            if cds_parameter[
                                0] == self.TRANSLATION_CDS_KEYWORD:
                                translation = [cds_parameter[1]]
                                while not current_line.endswith('"'):
                                    current_line = input_gbk_file.readline().strip()
                                    separate_line = current_line.split()
                                    translation.append(separate_line[0])
                                translation_line = ''.join(translation)
                                self.get_last_dict_value(parsed_data)[2][-1][2] \
                                    = translation_line.replace('"', '')
                                break
                            current_line = input_gbk_file.readline()
                            separate_line = current_line.split()

                current_line = input_gbk_file.readline().strip()
                while current_line != self.END_OF_PART_GBK_KEYWORD:
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
                                    left_translation \
                                        .append(
                                        parsed_data[current_key][2][k][2])
                                for j in range(i + 1, i + n_after + 1):
                                    if i >= len(parsed_data[current_key][2]):
                                        break
                                    right_translation \
                                        .append(
                                        parsed_data[current_key][2][j][2])
                                left_join = ''.join(left_translation[::-1])
                                right_join = ''.join(right_translation)
                                merged_translation = left_join + right_join
                                output_fasta_file.write(current_key + '\n' \
                                                        + merged_translation
                                                        + '\n')