from typing import Tuple


def fastq_filter(seqs:dict, gc_bounds:int | float | Tuple = (20,80),
                 length_bounds:int |float
                 |Tuple = (0, 2**32), quality_threshold:int = 0) -> dict:
    if seqs is None:
        raise ValueError('Your fastq_files are None')
    elif len(seqs) == 0:
        raise ValueError('Your fastq_files are empty')

    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int) or isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)

    seqs_filtered = {}
    for fastq_name in seqs:
        if int(seqs[fastq_name][-1]) >= quality_threshold and \
                length_bounds[0] <= len(seqs[fastq_name][0]) <= \
                length_bounds[1]:
            gc_content = ((seqs[fastq_name][0].count('G')
                          + seqs[fastq_name][0].count('C'))
                          / len(seqs[fastq_name][0]) * 100)
            if gc_bounds[1] >= gc_content >= gc_bounds[0]:
                seqs_filtered[fastq_name] = seqs[fastq_name]
    return seqs_filtered

