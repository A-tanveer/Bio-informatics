

def hamming_distance(a, b):
    """returns hamming distance between two strings."""
    if not len(a) == len(b):
        raise ValueError('Strings must be same in size to compute hamming distance')
    return sum(i != j for i, j in zip(a, b))


def ngram(string_to_gram, n):
    """returns list of ngrams of a strings"""
    n_gram = []
    for i in range(len(string_to_gram) - n + 1):
        n_gram.append(string_to_gram[i:i + n])
    return n_gram


def all_possible_kmers(k):
    """all kmers of DNA of length k"""
    from itertools import product
    alpha = ['A', 'C', 'G', 'T']
    return [''.join(i) for i in product(alpha, repeat=k)]


def d_neighbourhood(pattern, max_mismatch):
    """d neighbourhood
    :returns all possible kmers of the length of given pattern
    with maximum hamming distance of max_mismatch"""
    k_gram = all_possible_kmers(len(pattern))
    return [i for i in k_gram if hamming_distance(i, pattern) <= max_mismatch]


def median(list_of_dna, k):
    """takes a list of DNA of equal length and a number k.
    a list of median strings are returned. """
    kmers = all_possible_kmers(k)
    count_distance = []

    # checking hamming distance of every kmer with every gram of kgram of every DNA string
    # minimum hamming distance of a gram of each string for each kmer is stored.
    for x in kmers:
        count = 0
        for dna in list_of_dna:
            gram = ngram(dna, k)
            distance = []
            for each_gram in gram:
                distance.append(hamming_distance(x, each_gram))
            count += min(distance)
        count_distance.append(count)
    # then the kmers that have minimum hamming distance with all the strings are returned.
    m = min(count_distance)
    return [i for i, j in zip(kmers, count_distance) if j == m]


def random_dna(k):
    """Generate random DNA string of length k"""
    import random
    sign = 'ATGC'
    DNA = ''
    for i in range(k):
        DNA += sign[random.randint(0, 3)]
    return DNA


def random_rna(k):
    """Generate random RNA string of length k"""
    import random
    sign = 'AUGC'
    RNA = ''
    for i in range(k):
        RNA += sign[random.randint(0, 3)]
    return RNA


def num2dna(self, number, length=None):
    """converts number to a base 4 number. where digits are DNA letters.
    Adds A to the beginning of the dna if length of dna is less than length
    :returns the number converted into DNA sequence."""

    _alphabet = 'ACGT'
    _base = len(_alphabet)
    string = ''
    while number > 0:
        string = _alphabet[number % _base] + string
        number //= _base
    if length is not None:
        while len(string) < length:
            string = 'A' + string
    return string


def dna2num(self, string):
    """converts base 4 number(DNA sequence) to decimal number.
    :returns the DNA sequence converted into Decimal number"""

    _alphabet = 'ACGT'
    _base = len(_alphabet)
    number = 0
    for char in string:
        number = number * _base + _alphabet.index(char)
    return number
