from static import d_neighbourhood, all_possible_kmers, hamming_distance, ngram, median, random_dna


class BioInformatics:
    """@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    An approach to build a python library for bio informatics that covers the problems in Rosalind.info
    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"""

    def __init__(self, dna=''):
        """
        :param dna: initialize the DNA sequence
        """
        self.DNA = dna

    def implanted_motifs(self, dna_list, k, d):
        """takes a list of DNAs, length of kmers k, and maximum number of mismatch d
        :returns a list of implanted motifs"""

        kmers = all_possible_kmers(k)
        result = []
        for each in kmers:
            val = True
            for dna in dna_list:
                if len(self.approximate_matched_pattern(each, d, dna)) < 1:
                    val = False
                    break
            if val:
                result.append(each)
        return result

    def profile_most_probable_kmer(self, profile, dna=None):
        """takes a profile matrix. validate it, removes 0 error.
        :returns profile most probable kmer."""

        if dna is None:
            dna = self.DNA

        if True in [0 in pro for pro in profile]:
            print('profile contains 0. changing it...')
            for i in range(4):
                for j in range(len(profile[i])):
                    profile[i][j] += 0.2

        kmers = ngram(dna.upper(), len(profile[0]))
        pro = []
        for i in range(len(kmers)):
            b = list(kmers[i])
            x = 1
            for j in range(len(b)):
                latter = b[j]
                if latter == 'A':
                    x *= profile[0][j]
                if latter == 'C':
                    x *= profile[1][j]
                if latter == 'G':
                    x *= profile[2][j]
                if latter == 'T':
                    x *= profile[3][j]
            pro.append(x)
        return kmers[pro.index(max(pro))]

    def count_k_mer(self, sub_str, type='count', dna=None):
        """:returns the total number of occurrences of a kmer(sub_str) in the DNA string if type is count.
        :returns a list of starting indexes of the occurrences if type is index. """

        if dna is None:
            dna = self.DNA

        # can easily be done by using ngram
        # alternative simple implementation
        # count, starting_indexes = 0, []
        # kmers = ngram(dna, len(sub_str))
        # for i, j in enumerate(kmers):
        #     if j == sub_str:
        #         count += 1
        #         starting_indexes.append(i)
        # if type == 'index':
        #     return starting_indexes
        # return count

        starting_indexes = []
        if sub_str == '':
            return 0
        str_len, sub_len, start, count = len(dna), len(sub_str), sub_str[0], 0
        for k in range(str_len - sub_len + 1):
            if dna[k] == sub_str[0]:
                if dna[k:k + sub_len] == sub_str:
                    count += 1
                    starting_indexes.append(k)
        if type == 'index':
            return starting_indexes
        return count

    def reverse_complement(self, dna=None):
        """:returns reverse complement of a DNA sequence"""

        if dna is None:
            dna = self.DNA
        a = 'ATGCatgc'
        b = 'TACGtacg'
        complement = dna.translate({ord(x): y for (x, y) in zip(a, b)})
        return complement[::-1]

    def most_frequent_k_mer(self, k):
        """:returns the list of most frequent kmers in a DNA sequence."""
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

    def approximate_matched_pattern(self, pattern, max_mismatch, dna=None):
        """:returns a list of starting indexes of approximate matched patterns"""

        if dna is None:
            dna = self.DNA
        k_gram = ngram(dna, len(pattern))
        return [i for i, j in enumerate(k_gram) if hamming_distance(j, pattern) <= max_mismatch]

    def most_frequent_k_mer_with_mismatch(self, k, max_mismatch):
        """:returns the list of most frequent kmers of length k in
        a DNA sequence allowing maximum mismatch of max_mismatch"""

        k_mers = list(set(ngram(self.DNA, k)))
        frequency_counts = [len(self.approximate_matched_pattern(gram, max_mismatch)) for gram in k_mers]
        m = max(frequency_counts)
        return [k_mers[i] for i in range(len(frequency_counts)) if frequency_counts[i] == m]

    def most_frequent_k_mer_with_mismatch_and_complements(self, k, max_mismatch):
        """:returns the list of most frequent kmers of length k in a DNA sequence allowing
        maximum mismatch of max_mismatch also including occurrences of there reverse compliments"""

        k_mers = list(set(ngram(self.DNA, k)))
        frequency_counts = [(len(self.approximate_matched_pattern(gram, max_mismatch)) +
                             len(self.approximate_matched_pattern(self.reverse_complement(
                                 dna=gram), max_mismatch))) for gram in k_mers]
        m = max(frequency_counts)
        return [k_mers[i] for i in range(len(frequency_counts)) if frequency_counts[i] == m]

    def frequency_array(self, k, dna=None):
        """:returns frequency array of kmers of length k in DNA sequence"""
        if dna is None:
            dna = self.DNA
        kmers = all_possible_kmers(k)
        return [self.count_k_mer(kmer) for kmer in kmers]

    def lt_clump(self, k, l, t, dna=None):
        """
        :param k: length of kmers
        :param l: length L
        :param t:minimum number of kmers in L
        :param dna: optional DNA sequence
        :return: list of kmers for which L-t clump is present in the dna sequence.
        """
        if dna is None:
            dna = self.DNA
        clumps = []
        k_mers = list(set(ngram(dna, k)))
        for each in k_mers:
            for i in range(len(dna) - l + 1):
                if dna[i] == each[0]:
                    if self.count_k_mer(each, dna=dna[i:i + l]) >= t:
                        clumps.append(each)
                        break  # if present no need to check further for a kmer
        return clumps

    def skew(self, dna=None):
        """:returns: Y values of skew graph for dna sequence"""
        if dna is None:
            dna = self.DNA
        value = 0
        skew = [0]
        for i in dna:
            if i == 'C' or i == 'c':
                value -= 1
            elif i == 'g' or i == 'G':
                value += 1
            skew.append(value)
        return skew

    def minimum_skew(self, dna=None):
        """:returns X values of skew graph where Y is minimum"""
        if dna is None:
            dna = self.DNA
        skew = self.skew(dna)
        min_skew = []
        for i, j in enumerate(skew):
            if j == min(skew):
                min_skew.append(i)
        return min_skew
