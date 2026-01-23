## Programming Assignment #3

## do not remove the extra print statements
print('Question 1')
#################################################################
# QUESTION 1                                                    #
# Write a bit of Biopython code that will find the number of    # 
# occurrences and locations of ACYT within the following        # 
# sequence, where Y = C or T. Where are the matches located?    #
# Please note that you must only search the forward strand of   #
# the  sequence, which is provided on two lines here            # 
# Your code should read this sequence in as a fasta file that   #
# you create and upload with your solution. Print the positions #
# of the matches using 0-based positions, one per line          #
#                                                               #
# >sequence                                                     #
# AGCGATCTAGCATACTTATACGCGCGCAGCTATCGATCACTTGTGCTAGTAAAGTGCGC   #
# GCCGCATTAAAGTGCTAGCTAGCTACTTAGCTAGCTAGTCG                     #
#################################################################

##add code for question 1 here
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils
pattern = Seq("ACYT")

fastaFile = "seq_q1.fasta"
sequence_parse = SeqIO.parse(fastaFile, 'fasta')
for record in sequence_parse:
    print(record.seq)
    results = SeqUtils.nt_search(str(record.seq), pattern)
    print(results)

## do not remove the extra print statements
print('Question 2')
#################################################################
# QUESTION 2                                                    #
# Write a bit of Biopython code to transcribe the following DNA # 
# sequence and translate the resulting mRNA sequence. Your code #
# should print two lines of output, the transcribed RNA and     #
# the translated protein sequence                               #
#                                                               #
# >CDS                                                          #
# ATGCTGCGGCGAGCTCCGTCCAAAAGAAAATGGGGTTTGGTGTAAATCTGGGGGTGTAATG #
# TTATCATATAAATAAAAGAAAATGTAAAAAACAAAAACAAAAACAAAAAAGCC         #
#                                                               #
#################################################################

##add code for question 2 here
from Bio.Seq import Seq
CDS = Seq("ATGCTGCGGCGAGCTCCGTCCAAAAGAAAATGGGGTTTGGTGTAAATCTGGGGGTGTAATGTTATCATATAAATAAAAGAAAATGTAAAAAACAAAAACAAAAACAAAAAAGCC")
mRNA = CDS.transcribe()
print(mRNA)
protein = mRNA.translate()
print(protein)

## do not remove the extra print statements
print('Question 3')

#################################################################
# QUESTION 3                                                    #
# Write some Python code that will create and print out a       # 
# palindromic DNA sequence that is 200nt long. You can          # 
# start with random sequences that you type into the code       # 
# or that your code generates, but don't bother trying to type  #
# out a 200-nt palindromic DNA sequence by hand! Get creative   #
# how you  you create the sequence. Make sure you know what a   # 
# DNA palindrome is!                                            #                                   
#################################################################

##add code for question 3 here
from Bio.Seq import Seq
import random

nucleotides = "ATGC"
sequence = "".join(random.choices(nucleotides, k = 200))
x = Seq(sequence)
palindrome = x + x.reverse_complement()

if palindrome == palindrome.reverse_complement():
   print("Yes, it is a palindrome")
   print(palindrome) 




## do not remove the extra print statements
print('Question 4')
#################################################################
# QUESTION 4                                                    #
# For the following DNA sequences:                              #
#                                                               #
# GTAA                                                          #
# ACGT                                                          #
# GCGT                                                          #
# ACAT                                                          #
# GCGA                                                          #
#                                                               #
# Build a motif and weight matrix out of these sequences using  #
# biopython (1pt).  Use the resulting weight matrix to search   # 
# for matches in both strands of the following sequence.        #
#                                                               #
# >promoter_X1                                                  #
# ATGGTTGGCGTCACGAGCTTAGTCAGCCTAGCTAGCATAACGCGCTAGCAGCTAGCTACGC #
#                                                               #
# Your Biopython code should print out the position of where    #
# all the matches occur and their score. In other words, the    #
# output should be one line per match, with the position and    #
# score separated by a space. In other words:                   #
#                                                               #
# i1 S1                                                         #
# i2 S2                                                         #
# ...etc                                                        #
#                                                               #
# In addition, you should print out one more line with the      # 
# position and score of the highest scoring instance, so one    # 
# more line of output:                                          #
#                                                               #
# i_max S_max                                                   #
#                                                               #
# So your output should have one line per each instance/match   #
# and an additional line for the highest scoring match listed   #
# again on the last line of the output. If there are "ties" or  #
# equal scoring max positions, then print one of them in this   #
# format                                                        #  
#################################################################

##add code for question 4 here
from Bio.Seq import Seq
from Bio import motifs
sites = [Seq('GTAA'),Seq('ACGT'),Seq('GCGT'),Seq('ACAT'),Seq('GCGA')]
motif = motifs.create(sites)
print(motif.counts)
promoter_x1 = Seq('ATGGTTGGCGTCACGAGCTTAGTCAGCCTAGCTAGCATAACGCGCTAGCAGCTAGCTACGC')
k = len(motif.consensus)
motif.pseudocounts = 1
W = motif.pssm
i_max = None
S_max = float("-inf")

for i,S in W.search(promoter_x1):
    if i<0:
        pos = len(promoter_x1)+i
        print(pos,S)
    else:
        pos = i
        print(pos,S)
    
    if S > S_max:
        i_max = pos
        S_max = S
        
    

# Print best hit 
print("i_max S_max:", i_max, S_max)  

