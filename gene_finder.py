# -*- coding: utf-8 -*-
"""
Gene Coder

@author: Aurora Bunten

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
Salmonella = load_seq("./data/X73525.fa")


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
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
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse = ''
    for nucleotide in dna:
        complement = get_complement(nucleotide)
        reverse = reverse + complement
    return reverse[::-1]


get_reverse_complement('ATGCGT')


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
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
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
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
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
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
    another = find_all_ORFs(dna)
    rev = get_reverse_complement(dna)
    rev_another = find_all_ORFs(rev)
    both_rev = another + rev_another
    return(both_rev)


find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    total = find_all_ORFs_both_strands(dna)
    return(max(total, key=len))


longest_ORF("ATGCGAATGTAGCATCAAA")


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    long = []
    i = 0
    while i < num_trials:
        dna = shuffle_string(dna)
        long.append(longest_ORF(dna))
        i = i + 1
    return(max(long, key=len))

print(len(longest_ORF_noncoding(Salmonella, 5)))


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
    chain = []
    y = 0
    for y in range(0, len(dna), 3):
        codons = dna[y: (y+3)]
        new_aa = aa_table[codons]
        chain.append(new_aa)
    return chain


coding_strand_to_AA("ATGCGA")


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    amino = []
    threshold = longest_ORF_noncoding(dna, 1500)
    openORF = find_all_ORFs_both_strands(dna)
    for i in openORF:
        if len(i) < len(threshold):
            aa = coding_strand_to_AA(i)
            amino.append(aa)
    return aa

gene_finder(Salmonella)

if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
