# -*- coding: utf-8 -*-
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable

# Define the codon table for the organism
codon_table = CodonTable.unambiguous_dna_by_id[1]

# Define the gene regions and their corresponding start and end positions
gene_regions = {
    'gene1': {'start': 100, 'end': 500},
    'gene2': {'start': 600, 'end': 900},
    'gene3': {'start': 1100, 'end': 1500},
    'noncoding_region': {'start': 2000, 'end': 2500}
}

# Define a function to translate the nucleotide sequence to amino acids
def translate_sequence(sequence):
    """
    Translates a nucleotide sequence to an amino acid sequence
    using the standard genetic code
    """
    return Seq(sequence, generic_dna).translate(table=codon_table)

# Create a dictionary to map nucleotide positions to amino acid positions
nuc_to_aa = {}
for region, positions in gene_regions.items():
    if region == 'noncoding_region':
        # Skip noncoding regions
        continue
    start = positions['start']
    end = positions['end']
    nuc_seq = sequence[start-1:end]
    aa_seq = translate_sequence(nuc_seq)
    for i in range(len(nuc_seq)):
        nuc_pos = start + i
        aa_pos = i//3 + 1
        nuc_to_aa[nuc_pos] = {'region': region, 'aa_pos': aa_pos}

# Write the nuc_to_aa dictionary to the additional index file
with open('nuc_to_aa_index.txt', 'w') as f:
    for nuc_pos, info in nuc_to_aa.items():
        f.write('{}\t{}\t{}\t{}\n'.format(nuc_pos, info['region'], info['aa_pos'], info['aa_seq']))
