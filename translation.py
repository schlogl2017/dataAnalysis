#!usr/bin/env python


import codons
from codons import codon_table_std_RNA


def get_translation(sequence):
    codon = codons.get_codons(sequence, 'RNA')
    print(codon)


seq = 'ATGACAAGACTTCCTCTACTCAAACGACCTCGAAGAAATCGGAAAAGTGCTGCTATTCGATCTATGATTCGAGAAACTAA\
TATGGTTTCAAGTGATCTAATTTGGCCTATTTTTCTTAAAGAAGGTTCAGGAATTCGAGAAGAAATACCCAGTATGCCTG\
GAGTATACAGATGGAGTTTAGATACCATATCTAGAGAATTAGAAAGACTCTGTCTTATAGGGCTAAAAGCGGTTATCCTT\
TTCCCTGTTATAGAGGATCAGAAAAAAGATCAATTCGGAGCATATGCATCGCATCCGTATAACATTGTTTGTAGAGGAAT\
TCAGGAAATCAAAAAATCTTTTCCCCAACTGTGTGTGATTAGTGATATAGCTTTAGATCCTTTTACGACTAGTGGTCACG\
ATGGCATTTTTTACAACAATGAAGTGTTAAATGACGAAAGCGTTCGCGTA'