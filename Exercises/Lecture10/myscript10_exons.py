#!/usr/bin/python3

my_dna ="ATCGATCGATCGATCGACTGACTAGTCATAGCTATGCATGTAGCTACTCGATCGATCGATCGATCGATCGATCGATCGATCGATCATGCTATCATCGATCGATATCGATGCATCGACTACTAT"

# Remember that we start counting from zero, and that positions are inclusive at the start and exclusive at the end
exon1 = my_dna[0:63]
exon2 = my_dna[90:]

# Now we can just join the two exons and print them out
print("### Exons joined ###\n" + exon1 + exon2)

# To calculate the coding percentage, we just
# have to take the length of the exons,
# divide by the length of the sequence, and multiply by 100
coding_length = len(exon1 + exon2)
total_length = len(my_dna)
print("### Coding percentage (rounded) ###\n" + str(int((coding_length / total_length) * 100)))

# To print out the upper/lower case version,
# we need to take the middle bit
# i.e. the intron and convert it to lower case,
# then join the bits back together
# We already specified the exons
exon1 = my_dna[0:63]
exon2 = my_dna[90:]
intron = my_dna[63:90]
print("### Exon intron exon ###\n" + exon1 + intron.lower() + exon2)
