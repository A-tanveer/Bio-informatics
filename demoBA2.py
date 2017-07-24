from static import d_neighbourhood, all_possible_kmers, hamming_distance, ngram, median
import bio_informatics

x = bio_informatics.BioInformatics('')

# BA2A
# http://rosalind.info/problems/ba2a/
dnaStr = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
print(x.implanted_motifs(dnaStr, 3, 1))

# BA2B
# http://rosalind.info/problems/ba2b/
dnaStr = ['TGATGATAACGTGACGGGACTCAGCGGCGATGAAGGATGAGT',
          'CAGCGACAGACAATTTCAATAATATCCGCGGTAAGCGGCGTA',
          'TGCAGAGGTTGGTAACGCCGGCGACTCGGAGAGCTTTTCGCT',
          'TTTGTCATGAACTCAGATACCATAGAGCACCGGCGAGACTCA',
          'ACTGGGACTTCACATTAGGTTGAACCGCGAGCCAGGTGGGTG',
          'TTGCGGACGGGATACTCAATAACTAAGGTAGTTCAGCTGCGA',
          'TGGGAGGACACACATTTTCTTACCTCTTCCCAGCGAGATGGC',
          'GAAAAAACCTATAAAGTCCACTCTTTGCGGCGGCGAGCCATA',
          'CCACGTCCGTTACTCCGTCGCCGTCAGCGATAATGGGATGAG',
          'CCAAAGCTGCGAAATAACCATACTCTGCTCAGGAGCCCGATG']
print(median(dnaStr, 6))

# BA2C
# http://rosalind.info/problems/ba2c/
x.DNA = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
profile = [
    [0.2, 0.2, 0.3, 0.2, 0.3],
    [0.4, 0.3, 0.1, 0.5, 0.1],
    [0.3, 0.3, 0.5, 0.2, 0.4],
    [0.1, 0.2, 0.1, 0.1, 0.2]
]
print(x.profile_most_probable_kmer(profile))

# BA2D
# http://rosalind.info/problems/ba2d/

