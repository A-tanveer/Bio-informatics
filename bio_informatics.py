class BioInformatics:

    DNA = ''

    def __init__(self, DNA):
        self.DNA = DNA

    def ngram(self, string_to_gram, n):
        n_gram = []
        for i in range(len(string_to_gram) - n + 1):
            n_gram.append(self.DNA[i:i + n])
        return n_gram

    def k_mer(self, sub_str, type='count'):
        starting_indexes = []
        if sub_str == '':
            return 0
        str_len, sub_len, start, count = len(self.DNA), len(sub_str), sub_str[0], 0
        for k in range(str_len - sub_len + 1):
            if self.DNA[k] == sub_str[0]:
                if self.DNA[k:k+sub_len] == sub_str:
                    count += 1
                    starting_indexes.append(k)
        if type == 'index':
            return starting_indexes
        return count

    def most_frequent_k_mer(self, k):
        k_gram = self.ngram(self.DNA, k)
        set_gram = list(set(k_gram))
        set_count = []
        for a in set_gram:
            set_count.append(k_gram.count(a))
        m = max(set_count)
        result_list = []
        for a in range(len(set_count)):
            if set_count[a] == m:
                result_list.append(set_gram[a])
        return result_list

    def reverse_complement(self):
        a = 'ATGCatgc'
        b = 'TACGtacg'
        complement = self.DNA.translate({ord(x): y for (x, y) in zip(a, b)})
        return complement[::-1]

    def LT_clump(self, k, L, t):
        clumps = []
        k_mers = list(set(self.ngram(self.DNA, k)))
        for each in k_mers:
            for i in range(len(self.DNA) - L + 1):
                if self.DNA[i] == each[0]:
                    object1.DNA = self.DNA[i:i+L]
                    if object1.k_mer(each) >= t:
                        clumps.append(each)
                        break
        return clumps


object1 = BioInformatics('')
