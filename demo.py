import bio_informatics

# demo
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