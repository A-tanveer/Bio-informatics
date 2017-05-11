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


class BioInformatics:
    DNA = ''

    def __init__(self, DNA):
        self.DNA = DNA

    def k_mer(self, sub_str, type='count'):
        # can easily be done by using ngram
        starting_indexes = []
        if sub_str == '':
            return 0
        str_len, sub_len, start, count = len(self.DNA), len(sub_str), sub_str[0], 0
        for k in range(str_len - sub_len + 1):
            if self.DNA[k] == sub_str[0]:
                if self.DNA[k:k + sub_len] == sub_str:
                    count += 1
                    starting_indexes.append(k)
        if type == 'index':
            return starting_indexes
        return count

    def reverse_complement(self, dna=None):
        if dna is None:
            dna = self.DNA
        a = 'ATGCatgc'
        b = 'TACGtacg'
        complement = dna.translate({ord(x): y for (x, y) in zip(a, b)})
        return complement[::-1]

    def most_frequent_k_mer(self, k):
        result_list = []
        k_gram = ngram(self.DNA, k)
        set_gram = list(set(k_gram))
        set_count = []
        for a in set_gram:
            set_count.append(k_gram.count(a))
        m = max(set_count)
        for a in range(len(set_count)):
            if set_count[a] == m:
                result_list.append(set_gram[a])
        return result_list

    def approximate_matched_pattern(self, pattern, max_mismatch):
        k_gram = ngram(self.DNA, len(pattern))
        return [i for i in range(len(k_gram)) if hamming_distance(k_gram[i], pattern) <= max_mismatch]

    def most_frequent_k_mer_with_mismatch(self, k, max_mismatch):
        k_mers = list(set(ngram(self.DNA, k)))
        frequency_counts = [len(self.approximate_matched_pattern(gram, max_mismatch)) for gram in k_mers]
        m = max(frequency_counts)
        return [k_mers[i] for i in range(len(frequency_counts)) if frequency_counts[i] == m]

    def most_frequent_k_mer_with_mismatch_and_complements(self, k, max_mismatch):
        k_mers = list(set(ngram(self.DNA, k)))
        frequency_counts = [(len(self.approximate_matched_pattern(gram, max_mismatch)) +
                            len(self.approximate_matched_pattern(self.reverse_complement(dna=gram), max_mismatch)))
                            for gram in k_mers]
        m = max(frequency_counts)
        return [k_mers[i] for i in range(len(frequency_counts)) if frequency_counts[i] == m]

    def frequency_array(self, k, dna=None):
        if dna is None:
            dna = self.DNA
        kmers = all_possible_kmers(k)
        return [self.k_mer(kmer) for kmer in kmers]

    def num2dna(self, number, length=None):
        _alphabet = 'ACGT'
        _base = len(_alphabet)
        string = ''
        while number > 0:
            string = _alphabet[number % _base] + string
            number //= _base
        if not length is None:
            while len(string) < length:
                string = 'A' + string
        return string

    def dna2num(self, string):
        _alphabet = 'ACGT'
        _base = len(_alphabet)
        number = 0
        for char in string:
            number = number * _base + _alphabet.index(char)
        return number

    def LT_clump(self, k, L, t):
        clumps = []
        k_mers = list(set(ngram(self.DNA, k)))
        for each in k_mers:
            for i in range(len(self.DNA) - L + 1):
                if self.DNA[i] == each[0]:
                    object1.DNA = self.DNA[i:i + L]
                    if object1.k_mer(each) >= t:
                        clumps.append(each)
                        break
        return clumps

    def skew(self):
        temp = 0
        skew = [0]
        for i in self.DNA:
            if i == 'C' or i == 'c':
                temp -= 1
            elif i == 'g' or i == 'G':
                temp += 1
            skew.append(temp)
        return skew

    def minimum_skew(self):
        skew = self.skew()
        m = min(skew)
        min_skew = []
        for i in range(len(skew)):
            if skew[i] == m:
                min_skew.append(i)
        return min_skew


object1 = BioInformatics('')
