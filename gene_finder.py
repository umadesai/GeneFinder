# -*- coding: utf-8 -*-
"""
This is a gene finding program that can accurately determine regions of the Salmonella bacterium's DNA that code for proteins.

@author: Uma Desai

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
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

        Added doctests: Tests function for every possible nucleotide (this is not a strictly necessary doctest addition)
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    # TODO: implement this
    if nucleotide == "A":
    	return "T"
    elif nucleotide == "C":
    	return "G"
    elif nucleotide == "T":
    	return "A"
    elif nucleotide == "G":
    	return "C"


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

        These test that the function returns the complement of the DNA sequence backwards. 
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("CTGCCCGCGGG")
    'CCCGCGGGCAG'
    >>> get_reverse_complement("ATGCATGAATGTAGA")
    'TCTACATTCATGCAT'
    """
    # TODO: implement this
    dna_reversed = ''
    for i in dna:
    	dna_reversed += get_complement(i)
    return dna_reversed[::-1]

def rest_of_ORF(dna):
	""" Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

        Added doctests: Tests function to see if it returns the whole string when there is no in frame stop codon.
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGCATATATATATATATATATAT")
    'ATGCATATATATATATATATATAT'
    """
	have_found = False
	for i in range(3,len(dna)+1,3):
		if "TAG" == dna[i-3:i] or "TAA" == dna[i-3:i] or "TGA" == dna[i-3:i]:
			have_found = True
			break
	if have_found: 
		return dna[:i-3]
	return dna



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        These doctests test that nested ORFs are not included in the returned list of ORFs and that the sequence follows the multiples of three format
        (in the last test the stop codon is not in the correct multiple of three and thus does not end the frame.)
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGGACTGA")
    ['ATGGAC']
    >>> find_all_ORFs_oneframe("ATGGACCTGA")
    ['ATGGACCTGA']
    """
    non_nested_ORF = []
    i = 3
    while i < len(dna):
    	if dna[i-3:i] == "ATG":
    		rest_ORF = rest_of_ORF(dna[i-3:])
    		non_nested_ORF.append(rest_ORF)
    		i += len(rest_ORF)
    	i += 3
    return non_nested_ORF


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

        These test that if an ORF occurs entirely within another ORF it is not returned in the list that is outputed.
        
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC', 'ATGAATGTAGATAGATGTGCCC', 'ATG']
    >>> find_all_ORFs('ATGGACTGA')
    ['ATGGAC']
    """

    return find_all_ORFs_oneframe(dna[0:]) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

      	These test that the function checks for ORFs on both strands.

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC', 'ATGAATGTAGATAGATGTGCCC', 'ATG', 'ATGCAT']
    >>> find_all_ORFs_both_strands("ATGCGA")
    ['ATGCGA']
    """
    return find_all_ORFs(dna)+find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        These tests check that the function returns the longest ORF (the string of max length) of both strands of the dna sequence.

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGAATGTAG")
    'ATGAATGTAG'
    >>> longest_ORF("ATGCTACATTCGCAT")
    'ATGCTACATTCGCAT'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    return max(all_ORFs, key = len)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        There are no unit tests for this function because the shuffle_string function that is called gives a randomly ordered sequence of dna that can't be "expected" by a doctest.
        However, it is possible to check the function's output in the terminal and confirm that the results are reasonable.

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """

    long_ORFs = []
    for i in range(num_trials):
    	shuffled_dna = shuffle_string(dna)
    	ORFs_shuffled_dna = find_all_ORFs_both_strands(shuffled_dna)
    	long_ORFs.append(ORFs_shuffled_dna)
    return len(max(long_ORFs, key = len))


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        These doctests check that the function correctly computes the protein encoded by the dna sequence and concatenates and returns them in a string.

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("ATGCTTTTA")
        'MLL'
    """
    amino_acid_sequence = ''
    for i in range(3,len(dna)+1,3):
    	amino_acid = aa_table[dna[i-3:i]]
    	amino_acid_sequence += amino_acid
    return amino_acid_sequence

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

    """
    threshold = longest_ORF_noncoding(dna, 1500)
    ORFs = find_all_ORFs_both_strands(dna)
    sequences = []
    for i in range(len(ORFs)):
    	if len(ORFs[i]) > threshold:
    		sequence = coding_strand_to_AA(ORFs[i])
    		sequences.append(sequence)
    	else:
    		pass
    return sequences



if __name__ == "__main__":
	import doctest
    # doctest.run_docstring_examples(get_complement, globals())()
	# doctest.testmod()
	
	from load import load_seq
	dna = load_seq("./data/X73525.fa")
	print(gene_finder(dna))
	