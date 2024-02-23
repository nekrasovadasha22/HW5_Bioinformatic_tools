def phread33_converter(seqs: str) -> float | int:
    phread_score = 0
    for nucl_quality in seqs:
        phread_score += ord(nucl_quality) - 33
    return phread_score / len(seqs)


def convert_fastq_to_dict(input_path:str) -> dict:
    seqs = {}
    fastq_file_lines = []
    with open(input_path) as fastq_file:
        fastq_file_lines = fastq_file.readlines()
    for i in range(0, len(fastq_file_lines), 4):
        seqs[fastq_file_lines[i].strip()] = (fastq_file_lines[i + 1].strip(),
                                     fastq_file_lines[i + 3].strip())
    return seqs
