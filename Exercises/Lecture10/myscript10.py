#!/usr/bin/python3

# First try at calculating the A+T content in a DNA sequence
my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
# length = len(my_dna)
# a_count = my_dna.count('A')
# t_count = my_dna.count('T')
# at_content = (a_count + t_count) / length
# print("A+T content is " + str(at_content))

# The shorter version of the above
my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
at_content = (my_dna.count('A') + my_dna.count('T')) / len(my_dna)
print("A+T content is",str(int(100*at_content)),"percent")

# Complementing DNA
# replacement1 = my_dna.replace('A', 't')
# replacement2 = replacement1.replace('T', 'a')
# replacement3 = replacement2.replace('C', 'g')
# replacement4 = replacement3.replace('G', 'c')
# print("The complement sequence is " + str(replacement4.upper()))

# Let's try shortening the script into a one-liner
print("The complement sequence is", my_dna.replace('A', 't').replace('T','a').replace('G','c').replace('C','g').upper())

res = "GAATTC"
print("Lengths of",res,"cleaved fragments are",my_dna.find(res) + 1,"and",len(my_dna) -(my_dna.find(res)+1),"bases")


