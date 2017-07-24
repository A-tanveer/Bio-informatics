

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


def min_coins_num(amount, coins):
    """ Find the Minimum Number of Coins Needed to Make Change """
    num_list = [0]*(amount+1)
    num_list[0] = 0
    for k in range(1, amount+1):
        cur_min = float('inf')
        for i in coins:
            if k >= i:
                cur_min = min(cur_min, num_list[k-i])
                print(cur_min)

        num_list[k] = cur_min + 1
        print(num_list)
        print()

    print(num_list[-1])
    return num_list[-1]




def longest_com_subseq(s1, s2):
    """
    :param s1:
    :param s2:
    :return: longest common substring
    """
    D = [[0]*(len(s2)+1) for _ in range(len(s1)+1)]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                D[i+1][j+1] = D[i][j]+1
            else:
                D[i+1][j+1] = max(D[i+1][j], D[i][j+1])
    sub_seq = ''
    i = len(s1)
    j = len(s2)
    while i != 0 and j != 0:
        if D[i][j] == D[i-1][j]:
            i -= 1
        elif D[i][j] == D[i][j-1]:
            j -= 1
        else:
            sub_seq = s1[i-1] + sub_seq
            i -= 1
            j -= 1
    return sub_seq


if __name__ == '__main__':
    # n, coins = read_data('in.txt')
    # min_coins_num(40, [1,5,10,20,25,50])
    print(longest_com_subseq('PLEASANTLY', 'MEANLY'))