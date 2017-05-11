import bio_informatics

x = bio_informatics.BioInformatics('ACAACTATGCATACTATCGGGAACTATCCT')

#BA1A http://rosalind.info/problems/ba1a/
print(x.k_mer('ACTAT')) # = 3 k-mers

# BA1B http://rosalind.info/problems/ba1b/
print(x.most_frequent_k_mer(4))

# BA1C http://rosalind.info/problems/ba1c/
print(x.reverse_complement())

# BA1D http://rosalind.info/problems/ba1d/
x.DNA = 'GATATATGCATATACTT'
print(x.k_mer('ATAT', type='index'))

# BA1E http://rosalind.info/problems/ba1e/
x.DNA = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
print(x.LT_clump(5, 75, 4))

# BA1F http://rosalind.info/problems/ba1f/
x.DNA = 'CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG'
print(x.minimum_skew())

# BA1G http://rosalind.info/problems/ba1g/
