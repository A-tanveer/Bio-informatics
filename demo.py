import bio_informatics

# demo
x = bio_informatics.BioInformatics('ACAACTATGCATACTATCGGGAACTATCCT')

#BA1A http://rosalind.info/problems/ba1a/
print(x.k_mer('ACTAT')) # = 3 k-mers

# BA1B http://rosalind.info/problems/ba1b/
print(x.most_frequent_k_mer(4))

