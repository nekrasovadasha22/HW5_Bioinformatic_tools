def phread33_converter(seqs: str) -> float | int:
    phread_score = 0
    for nucl_quality in seqs:
        phread_score += ord(nucl_quality) - 33
    return phread_score / len(seqs)
