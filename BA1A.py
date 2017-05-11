
# http://rosalind.info/problems/ba1a/


class BioInformatics:
    def __init__(self):
        pass

    def k_mer(self, string, sub_str):
        if sub_str == '':
            return 0

        str_len, sub_len, start, count = len(string), len(sub_str), sub_str[0], 0

        for k in range(str_len - sub_len + 1):
            if string[k] == sub_str[0]:
                if string[k:k+sub_len] == sub_str:
                    count += 1
        return count

# demo
x = BioInformatics()
print(x.k_mer('ACAACTATGCATACTATCGGGAACTATCCT', 'ACTAT')) # = 3 k-mers
