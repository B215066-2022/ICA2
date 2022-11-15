#!/usr/bin/python3

my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
length = len(my_dna)
a_count = my_dna.count('A')
t_count = my_dna.count('T')
at_content = (a_count + t_count) / length
print("A+T content is " + str(at_content))

# Program to take a DNA sequence and calculate A+T content
# Written by s123456 on 21 Oct 2022
# -------------------------------------------- #
# Input
my_dna = "ACTGATCGATTACGTATAGTATTTGCTATCATACATATATATCGATGCGTTCAT"
# Process
at_content = (my_dna.count('A') + my_dna.count('T')) / len(my_dna)
# Output
print("A+T content is",str(int(100*at_content)),"percent")


