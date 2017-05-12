from static import d_neighbourhood, all_possible_kmers, hamming_distance, ngram
import bio_informatics

x = bio_informatics.BioInformatics('')

# BA2A
# http://rosalind.info/problems/ba2a/
dnaStr = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
print(x.implanted_motifs(dnaStr, 3, 1))
