class BioInformatics:

    DNA = ''

    def __init__(self, DNA):
        self.DNA = DNA
        pass

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

        k_gram = []

        for i in range(len(self.DNA) - k + 1):
            k_gram.append(self.DNA[i:i+k])

        set_gram = set(k_gram)
        set_gram = list(set_gram)
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

