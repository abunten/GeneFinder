# -*- coding: utf-8 -*-
"""
Gene Coder

@author: Aurora Bunten

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'T':
        return 'A'


get_complement('C')


def get_reverse_complement(dna):
    reverse = ''
    for nucleotide in dna:
        complement = get_complement(nucleotide)
        reverse = reverse + complement
    return reverse[::-1]


get_reverse_complement('ATGCGT')


def rest_of_ORF(dna):
    stop1 = 'TGA'
    stop2 = 'TAG'
    stop3 = 'TAA'
    y = 0
    for y in range(0, len(dna), 3):
        if dna[y: (y+3)] == stop1:
            rest = dna[:y]
            return rest
        if dna[y: (y+3)] == stop2:
            rest = dna[:y]
            return rest
        if dna[y: (y+3)] == stop3:
            rest = dna[:y]
            return rest
    return dna


rest_of_ORF("ATGAGATAGG")


def find_all_ORFs_oneframe(dna):
    i = 0
    ORF = []
    while (i < len(dna)):
        if dna[i: i+3] == 'ATG':
            thisORF = rest_of_ORF(dna[i:])
            ORF.append(thisORF)
            i = i + len(thisORF)
        else:
            i = i+3
    return ORF


find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")


def find_all_ORFs(dna):
    one = find_all_ORFs_oneframe(dna)
    two = find_all_ORFs_oneframe(dna[1:])
    three = find_all_ORFs_oneframe(dna[2:])
    total = []
    total.extend(one)
    total.extend(two)
    total.extend(three)
    return total


find_all_ORFs("ATGCATGAATGTAGATAGATGTGCCC")

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
