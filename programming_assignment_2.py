
## do not remove the extra print statements
print('Question 1')
#################################################################
# QUESTION 1                                                    #
# For the DNA sequence "TCATCTGCCCGGCAT", create a Seq object   #
# called DNA for this sequence, then make a second variable     #
# called RNA storing the corresponding RNA. Then print the      #
# reverse complement of the RNA. Finally, print the translation #
# of the reverse complement of the RNA. Please do not print any #
# additional text other than the sequences described above.     #
# Additonal text, such as "RNA:" will be misread by the autograder
# so only print the sequences requested 
#################################################################



##add code for question 1 here
from Bio.Seq import Seq
DNA = Seq("TCATCTGCCCGGCAT")
RNA = DNA.transcribe()
DNA_RC = DNA.reverse_complement()
RNA_RC = DNA_RC.transcribe()
print (RNA_RC)
protein_rc = RNA_RC.translate()
print(protein_rc)



## do not remove the extra print statements
print('Question 2')
#################################################################
# QUESTION 2                                                    #
# For the DNA sequence "ATCAGAGCTAGCTAGCAC", create a Seq       #
# object with this sequence, then print the substring that      #
# starts at the first "G" and is of length 3. Only print the    #
# subsequence/substring and nothing else in the section below:  #
#################################################################


##add code for question 2 here
from Bio.Seq import Seq
DNA = Seq("ATCAGAGCTAGCTAGCAC")
sub = str(DNA)
i = sub.index("G")
print(DNA[i:i+3])



## do not remove the extra print statements
print('Question 3')
#################################################################
# QUESTION 3                                                    #
# You can print/create the complementary sequence for any DNA   #
# sequence using the complement() method as part of biopython   #
# Seq object. The complementary sequence swaps each nucleotide  #
# with its complementary character (e.g. G for C, A for T and   #
# so on) For example,                                           #
#                                                               #
#>>> from Bio.Seq import Seq                                    #
#>>> DNA = Seq('GTCTGACGTA')                                    #
#>>> print(DNA.complement())                                    #
#CAGACTGCAT                                                     #
#                                                               #
# Write some code that will create the DNA sequence             #
# "GTCGTAGTGACTCA" as a Seq object, then print its complementary#
# sequence using Biopython methods associated with the Seq      #
# object, and without using the complement() method. Hint: you  #
# should think of reverse compelement as a two-step process     #
# (reverse and complement, I have a youtube video on this idea) #
# don't use the method shown in BB 345 where you create a comp  #
# dictionary, no need to. Don't overthink it! This is just a    #
# couple lines of code. Only print out the complementary DNA    #
#################################################################


##add code for question 3 here
from Bio.Seq import Seq
DNA = Seq("GTCGTAGTGACTCA")
rev_comp_DNA = DNA.reverse_complement()
comp_DNA = rev_comp_DNA[::-1]
print(comp_DNA)


## do not remove the extra print statements
print('Question 4')
#################################################################
# QUESTION 4                                                    #
# Write some code that will create the DNA sequence             # 
# "TTATCGTACGATGCTAGGCATC" as a Seq object, then print a        #
# subsequence that includes every third character starting at   #
# the first "A", and is of length 4. Only print the subsequence #
#################################################################


##add code for question 4 here
from Bio.Seq import Seq
DNA = Seq("TTATCGTACGATGCTAGGCATC")
sub = str(DNA)
i =sub.index("A")
sub_DNA = DNA[i::3]
final_sub_dna = sub_DNA[:4]
print(final_sub_dna)



## do not remove the extra print statements
print('Question 5')
#################################################################
# QUESTION 5                                                    #
# There is a FASTA file of fungal 18S RNAs included with this   #
# assignment called 18S_rRNAs.fasta. This file is referenced    #
# by a variable in the starter code below as fasta_file         #
# Copy this file to your path, and write a python script that   #
# will read the FASTA file, and for each of these sequences,    #
# print out the reverse complement of last 50 nucleotides. You  #
# should first extract the last 50 nucleotides, and then print  #
# the reverse complement of that subsequence.                   #
# One important point. If you test this on the servers, the path#
# to the file will be different than where it is on the Autograder
# The path below is for the Autograder, but you may want to 
# *temporarily* change the path for testing out on the server
# but change it back before submitting to the autograder. 
#################################################################


from Bio import SeqIO
fasta_file = '/autograder/source/18S_rRNAs.fasta'
##add code for question 5 here
sequences = SeqIO.parse(fasta_file, 'fasta')
for record in sequences:
    rec_id = record.id
    seq = record.seq
    final_50 = seq[-50:]
    rc_final_50 = final_50.reverse_complement()
    print(rec_id, rc_final_50)
