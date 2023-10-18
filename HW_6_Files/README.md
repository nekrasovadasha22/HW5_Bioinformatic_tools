# Functions for Handling Bioinformatics Files

This Python module provides two utility functions for handling bioinformatics file formats, specifically FASTA and GenBank (GBK) files. You can use these functions to streamline your data processing tasks.

## Function 1: convert_multiline_fasta_to_oneline

### Description
This function is used to convert a multi-line FASTA file into a one-line FASTA file.

### Usage
```python
from bio_files_processor import convert_multiline_fasta_to_oneline
```

# Example usage:

```python
convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
```
`input_fasta` (str): Path to the input multi-line FASTA file.

`input_fasta` = "input.fasta"
```python
'''>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
GGTGCTGTG'''
```
`output_fasta` (str, optional): Path to the output one-line FASTA file (default is generated if not provided
convert_multiline_fasta_to_oneline(input_fasta, output_fasta)
`output_fasta` = "output.fasta"
```python
'''>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+) 
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG'''
```

## Function 2: select_genes_from_gbk_to_fasta
### Description
This function is used to select specific genes from a GenBank (GBK) file and save them protein sequence in a FASTA file.
### Usage
```python
from bio_files_processor import select_genes_from_gbk_to_fasta
```
# Example usage:
`input_gbk` = "input.gbk"
`input_gbk` (str): Path to the input GenBank file.

`genes` (list): A list of gene names to select.

`genes` = ["GeneA", "GeneB"]

`output_fasta` (str, optional): Path to the output FASTA file (default is generated if not provided).

`output_fasta` = "output.fasta"

`n_before` (int, optional): Number of genes before the interesting gene (default is 1).

`n_after` (int, optional): Number of genes after the interesting gene (default is 1).
```python
select_genes_from_gbk_to_fasta(input_gbk, ['dtpD'])

#input_gbk_file
'''FEATURES             Location/Qualifiers
     source          1..2558431
                     /organism="Escherichia coli"
                     /mol_type="genomic DNA"
                     /strain="strain"
                     /db_xref="taxon:562"
     CDS             384..893
                     /locus_tag="IFLAKNEJ_00001"
                     /inference="ab initio prediction:Prodigal:002006"
                     /codon_start=1
                     /transl_table=11
                     /product="hypothetical protein"
                     /translation="MNQQRFDDSTLIRIFALHELHRLKEHGLTRGALLDYHSRYKLVF
                     LAHSQPEYRKLGPFVADIHQWQNLDDFYNQYYQRVIVLLSHPANPRDHTNVLMHVQGY
                     FRPHIDSTERQQLAALIDSYRRGEQPLLAPLMRIKHYMALYPDAWLSGQRYFELWPRV
                     INLRHSGVL"
     CDS             890..2308
                     /gene="phrB"
                     /locus_tag="IFLAKNEJ_00002"
                     /EC_number="4.1.99.3"
                     /inference="ab initio prediction:Prodigal:002006"
                     /inference="similar to AA sequence:UniProtKB:P00914"
                     /codon_start=1
                     /transl_table=11
                     /product="Deoxyribodipyrimidine photo-lyase"
                     /db_xref="COG:COG0415"
                     /translation="MITHLVWFRQDLRLHDNLALAAACRNSSARVLALYIATPRQWAT
                     HNMSPRQAELINAQLNGLQIALAEKGIPLLFREVDDFVASVEIVKQVCAENSVTHLFY
                     NYQYEVNERARDVEVERALRNVVCEGFDDSVILPPGAVMTGNHEMYKVFTPFKNAWLK
                     RLREGMPECVAAPKVRSSGSIEPAPSITLNYPRQSFDTAHFPVEEKAAIAQLRQFCQN
                     GAGEYEQQRDFPAVEGTSRLSASLATGGLSPRQCLHRLLAEQPQALDGGAGSVWLSEL
                     IWREFYRHLMTYYPSLCKHCPFIAWTDRVQWQSNPAHLQAWQKGKTGYPIVDAAMRQL
                     NSTGWMHNRLRMITASFLVKDLLIDWREGERYFMSQLIDGDLAANNGGWQWAASTGTD
                     AAPYFRIFNPITQGEKFDREGEFIRRWLPELRDVPGKAVHEPWKWAQKAGVKLDYPQP
                     IVDHKEARLRTLAAYEEARKGA"
     CDS             complement(2350..3831)
                     /gene="dtpD"
                     /locus_tag="IFLAKNEJ_00003"
                     /inference="ab initio prediction:Prodigal:002006"
                     /inference="similar to AA sequence:UniProtKB:P75742"
                     /codon_start=1
                     /transl_table=11
                     /product="Dipeptide permease D"
                     /db_xref="COG:COG3104"
                     /translation="MNKHASQPRAIYYVVALQIWEYFSFYGMRALLILYLTNQLKYND
                     THAYELFSAYCSLVYVTPILGGFLADKVLGNRMAVMLGALLMAIGHVVLGASEIHPSF
                     LYLSLAIIVCGYGLFKSNVSCLLGELYEPTDPRRDGGFSLMYAAGNVGSIIAPIACGY
                     AQEEYSWAMGFGLAAVGMIAGLVIFLCGNRHFTHTRGVNKKVLRATNFLLPNWGWLLV
                     LLVATPALITVLFWKEWSVYALIVATIIGLGVLAKIYRKAENQKQRKELGLIVTLTFF
                     SMLFWAFAQQGGSSISLYIDRFVNRDMFGYTVPTAMFQSINAFAVMLCGVFLAWVVKE
                     SVAGNRTVRIWGKFALGLGLMSAGFCILTLSARWSAMYGHSSLPLMVLGLAVMGFAEL
                     FIDPVAMSQITRIEIPGVTGVLTGIYMLLSGAIANYLAGVIADQTSQASFDASGAINY
                     SINAYIEVFDQITWGALACVGLVLMIWLYQALKFRNRALALES"'''
#output_gbk_file
'''
NODE_1_length_2558431_cov_75.1851642558431
MITHLVWFRQDLRLHDNLALAAACRNSSARVLALYIATPRQWATHNMSPRQAELINAQLNGLQIALAEKGIPLLFREVDDFVASVEIVKQVCAENSVTHLFYNYQYEVNERARDVEVERALRNVVCEGFDDSVILPPGAVMTGNHEMYKVFTPFKNAWLKRLREGMPECVAAPKVRSSGSIEPAPSITLNYPRQSFDTAHFPVEEKAAIAQLRQFCQNGAGEYEQQRDFPAVEGTSRLSASLATGGLSPRQCLHRLLAEQPQALDGGAGSVWLSELIWREFYRHLMTYYPSLCKHCPFIAWTDRVQWQSNPAHLQAWQKGKTGYPIVDAAMRQLNSTGWMHNRLRMITASFLVKDLLIDWREGERYFMSQLIDGDLAANNGGWQWAASTGTDAAPYFRIFNPITQGEKFDREGEFIRRWLPELRDVPGKAVHEPWKWAQKAGVKLDYPQPIVDHKEARLRTLAAYEEARKGAMKNTELEQLINEKLNSAAISDYAPNGLQVEGKETVQKIVTGVTASQALLDEAVRLGADAVIVHHGYFWKGESPVIRGMKRNRLKTLLANDINLYGWHLPLDAHPELGNNAQLAALLGITVMGEIEPLVPWGELTMPVPGLELASWIEARLGRKPLWCGDTGPEVVQRVAWCTGGGQSFIDSAARFGVDAFITGEVSEQTIHSAREQGLHFYAAGHHATERGGIRALSEWLNENTDLDVTFIDIPNPA
```