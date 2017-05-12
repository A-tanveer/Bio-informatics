def hamming_distance(a, b):
    if not len(a) == len(b):
        raise ValueError('Strings must be same in size to compute hamming distance')
    return sum(i != j for i, j in zip(a, b))


def ngram(string_to_gram, n):
    n_gram = []
    for i in range(len(string_to_gram) - n + 1):
        n_gram.append(string_to_gram[i:i + n])
    return n_gram


def all_possible_kmers(k):
    from itertools import product
    alpha = ['A', 'C', 'G', 'T']
    return [''.join(i) for i in product(alpha, repeat=k)]


def d_neighbourhood(pattern, max_mismatch):
    k_gram = all_possible_kmers(len(pattern))
    return [i for i in k_gram if hamming_distance(i, pattern) <= max_mismatch]


def median(dna, k):
    kmers = all_possible_kmers(k)
    count_distance = []
    for x in kmers:
        count = 0
        for str in dna:
            gram = ngram(str, k)
            distance = []
            for each in gram:
                distance.append(hamming_distance(x, each))
            count += min(distance)
        count_distance.append(count)
    m = min(count_distance)
    return [i for i,j in zip(kmers, count_distance) if j == m]
