# Utility Functions for Handling Bioinformatics Files

This Python module provides two utility functions for handling bioinformatics file formats, specifically FASTA and GenBank (GBK) files. You can use these functions to streamline your data processing tasks.

## Function 1: convert_multiline_fasta_to_oneline

### Description
This function is used to convert a multi-line FASTA file into a one-line FASTA file.

### Usage
`python
from biofileutils import convert_multiline_fasta_to_oneline

# Example usage:
input_fasta = "input.fasta"
output_fasta = "output.fasta"
convert_multiline_fasta_to_oneline(input_fasta, output_fasta)

`input_fasta` (str): Path to the input multi-line FASTA file. 


`output_fasta` (str, optional): Path to the output one-line FASTA file (default is generated if not provided
