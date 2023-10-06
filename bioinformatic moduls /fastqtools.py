from typing import Tuple


def phread33_converter(seqs):
    phread_score = 0
    for nucl_quality in seqs:
        phread_score = (ord(nucl_quality) - 33) / len(nucl_quality)
    return phread_score




