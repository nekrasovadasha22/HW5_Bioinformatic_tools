from typing import Tuple
import fastqtools

img_fasta_file = {'my_fasta':('ATGCUGCGCT', '33')}

fastqtools.fastq_filter(img_fasta_file)
print(fastqtools.fastq_filter(img_fasta_file))