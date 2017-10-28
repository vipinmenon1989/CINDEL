
This code is developed by A Vvipin Menon, and BIG lab retains all the rights concerned with code 

This code is built for Cpf1 based genome editing to identify highly active guide RNAs from a set of guide RNAs.

 single code for calcualting score for CRISPR-Cpf1 sgRNA,CINDEL score 
# The code  comprises of Free Energy ,Mononucleotide ,Dinucleotide and Independent nucleotide composition
# This code generates two results one for whole batch and other one for single sequence,option is given 
# The total length of the sequence should be  27 bp , whereas the target sequence is 23 bp and 4bp would be PAM , Please comply with PAM : TTTV (V = A,G or C) 

The basic example of running the code is 

a) if you want to run for ordered set of sequences and generate output in file with score 
	
	python CINDEL.py -a input file

b) If you want to generate just score for a sequence

	python CINDEL.py -b Sequence
c) If you want to find potential Cpf1 sequences within a set sequence
	
	python CINDEL.py -c guidefinder
 