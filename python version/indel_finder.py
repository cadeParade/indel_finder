from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np
import string
import re
from decimal import *

#epas

wt = "cactgttgttaggagggttcagagtagcaggatgaagttgctgttgtttattttggatgtgagccaagggcttgagagctagaacaagactagtatagtgtgcacacacactaacttgcattctaaaactcttgtgtttgtgctgtattgcagGCTACAATACTCCCACTGAAATGACAGATGCAGACAGACTCATGGACAGTTGGTATCTGAAGTCACTCGGTGGCTTTATTACAGTGGTAACATCAGATGGAGACATGATCTTCTTATCGGAGAACATCAACAAtagtaacgcacactgtatcaacacatgaatcga"

a6 = "NNNNNNNNNNNNGANNCTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTCACTGTTGTTAGGAGGGTTCAGAGTAGCAGGATGAAGTTGCTGTTGTTTATTTTGGATGTGAGCCAAGGGTTTGAGAGCTAGAACAAGACTAGTATAGTGTGCACACACACTAACATGCATTCTAAAACTCTTGTGTTTGTGCTGTATTGCAGGCTACGCTACTCCCACTGAAATGACAGACTCATCCATGACAGACTCATGGACAGTTGGTATCTGAAGTCACTCGGTGGCTTTATTACAGTGGTAACATCAGATGGAGACATGATCTTCTTATCGGAGAACATCAACAAATTCATGGGTCTCACTCAGGTGAGTAGTAACGCACACTGTATCAACACATGAATCGAAATCCCGCGGCCATGGCGGCCGGGAGCATGCGACGTCGGGCCCAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGACGCGCCCTGTANCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCNTANCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCNCNNNCNNNCNNNNCTTTNCCNNNCNNGCTNNNANCGGGGGCTNCCTTTNNGGGNTNCCNATTTANNNCTTTANNGCNNCNNNNNNNNNNNNN"



bind1 = "TACAATACTCCCACTGAA"
bind2 = "ACTCATGGACAGTTGGTA"
target = "ATGACAGATGCAGACAG"


def read_fasta(file_to_open):
	
	seq_file = open(file_to_open, "rU")
	sequence_list = []
	
	for record in SeqIO.parse(seq_file, "fasta", 
							  alphabet = generic_dna):
		sequence_list.append(record)
		
	seq_file.close()
	return sequence_list
	
def annotate_sequence(sequence, bpcutoff, npctcutoff):

	#assign well number to name field
	current_well_no = re.split("([A-H])(\\d\\d)", sequence.name)
	current_well_no = "".join(current_well_no[1:3])
	sequence.name = current_well_no

	#assign length and n count
	sequence.length = len(sequence)
	sequence.n_count = sequence.seq.count("N")
	
	#find percent "N" of sequence
	sequence.n_pct = (Decimal(sequence.n_count)/
					  Decimal(sequence.length)) * 100
	
	#screens for too few bps or too many N's
	if (sequence.n_pct > npctcutoff and 
		sequence.length < bpcutoff):
		sequence.initial_screen = ("sequence too short" + 
								   "and too many N's")
	elif sequence.n_pct > npctcutoff:
		sequence.initial_screen = "too many N's"
	elif sequence.length < bpcutoff:
		sequence.initial_screen = "sequence too short"
	else:
		sequence.initial_screen = "OK"	
	
	
# def annotate_sequences(sequence_list, bpcutoff, npctcutoff):
# 	
# 	sequence_dict = {}
# 	first_binding_site_indexes = []
# 	second_binding_site_indexes = []
# 	first_binding_site_indexes_rev = []
# 	second_binding_site_indexes_rev = []
# 		
# 	for i,sequence in enumerate(sequence_list):
	
		# #assign well number to name field
# 		current_well_no = re.split("([A-H])(\\d\\d)", sequence.name)
# 		current_well_no = "".join(current_well_no[1:3])
# 		sequence.name = current_well_no
# 		
# 		#assign length and n count
# 		sequence.length = len(sequence)
# 		sequence.n_count = sequence.seq.count("N")
# 		
# 		#find percent "N" of sequence
# 		sequence.n_pct = (Decimal(sequence.n_count)/
# 						  Decimal(sequence.length)) * 100
# 		
# 		#screens for too few bps or too many N's
# 		if (sequence.n_pct > npctcutoff and 
# 			sequence.length < bpcutoff):
# 			sequence.initial_screen = ("sequence too short" + 
# 									   "and too many N's")
# 		elif sequence.n_pct > npctcutoff:
# 			sequence.initial_screen = "too many N's"
# 		elif sequence.length < bpcutoff:
# 			sequence.initial_screen = "sequence too short"
# 		else:
# 			sequence.initial_screen = "OK"	
# 		
#		first_binding_site_indexes.append(
#				find_best_match_location(sequence, bind1))
#		first_binding_site_indexes_rev.append(
	# 				find_best_match_location(
# 					sequence.reverse_complement(), bind1))
# 		print sequence.name
# 		print first_binding_site_indexes[i][2], "forward"
# 		print first_binding_site_indexes_rev[i][2], "reverse"
# 		
# 	
	# print sequence_list[1].name, "name"
# 	print sequence_list[1].length, "length"
# 	print sequence_list[1].n_count, "n count"
# 	print sequence_list[1].n_pct, "n pct"
# 	print sequence_list[1].initial_screen, "initial screen"
# 	print sequence_list[1]
	
def decide_sequence_direction(sequence, fwd_tuple, rev_tuple):
	if fwd_tuple[2] < rev_tuple[2] and fwd_tuple[2] < 3:
		sequence.direction = "forward"
		sequence.bind1 = fwd_tuple
	elif fwd_tuple[2] > rev_tuple[2] and rev_tuple[2] < 3:
		sequence.direction = "reverse"
		sequence.bind1 = rev_tuple
	else:
		sequence.direction = "dunno"
	print sequence.name
	print sequence.direction
	print fwd_tuple[2]
	print rev_tuple[2]
	print sequence.n_pct
	print sequence.initial_screen
		


def find_best_match_location(sequence, search_string):
	
	score_list = []
	
	for i in range(0,len(sequence)-len(search_string)):
		current_string = sequence[i: i+len(search_string)]
		score = 0
		for i,base in enumerate(current_string):
			if base == search_string[i]:
				score = score
			else:
				score +=1
		score_list.append(score)
	
	min_score_index = score_list.index(min(score_list))
	
	#returns tuple of search strings starting and ending index
	# and match score
	return (min_score_index, 
			min_score_index + len(search_string), 
			min(score_list))
	

def main():
	
	#### SETUP ######
	file_to_open = "epas1b.seq"
	bp_cutoff_value = 275
	n_pct_cutoff_value = 32
	#sets how many digits for numbers defined as Decimal(#)
	getcontext().prec = 3
	
	#running	
	sequence_list = read_fasta(file_to_open)
	
	for i, sequence in enumerate(sequence_list):
		annotate_sequence(sequence, 
						  bp_cutoff_value, 
						  n_pct_cutoff_value)
		fwd  = find_best_match_location(sequence, bind1)
		rev = find_best_match_location(sequence.reverse_complement(),
									   bind1)
		decide_sequence_direction(sequence, fwd, rev)
		
	
	#annotate_sequences(sequence_list,
	#				   bp_cutoff_value, 
	#				   n_pct_cutoff_value)
					   
					   
	
main()