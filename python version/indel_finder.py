from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np
import string
import re
from decimal import *

#nrxn
bind1 = "ATCTTCAGC" 
bind2 = "GATGAGGTT"

#epas
# bind1 = "TACAATACTCCCACTGAA"
# bind2 = "ACTCATGGACAGTTGGTA"
# target = "ATGACAGATGCAGACAG"
#wt = "cactgttgttaggagggttcagagtagcaggatgaagttgctgttgtttattttggatgtgagccaagggcttgagagctagaacaagactagtatagtgtgcacacacactaacttgcattctaaaactcttgtgtttgtgctgtattgcagGCTACAATACTCCCACTGAAATGACAGATGCAGACAGACTCATGGACAGTTGGTATCTGAAGTCACTCGGTGGCTTTATTACAGTGGTAACATCAGATGGAGACATGATCTTCTTATCGGAGAACATCAACAAtagtaacgcacactgtatcaacacatgaatcga"


def read_fasta(file_to_open):
	
	seq_file = open(file_to_open, "rU")
	sequence_list = []
	
	for record in SeqIO.parse(seq_file, "fasta", 
							  alphabet = generic_dna):
		sequence_list.append(record)
		record.length = None 
		record.n_count = None 
		record.n_pct = None 
		record.initial_screen = "Omitted" 
		record.direction = None 
		record.bind1dir = None 
		record.bind2dir = None 
		record.distance_between_sites = None
		
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
		
def find_best_match_location(sequence, search_string):
	
	score_list = []
	#print sequence
	#print sequence.initial_screen
	#print len(sequence), "length of sequence"
	#print len(search_string), "length of search string"
	#if sequence.initial_screen == "OK":
	# Make list of scores for each window
	for i in range(0,len(sequence)-len(search_string)):
		current_string = sequence[i: i+len(search_string)]
		score = 0
		for i,base in enumerate(current_string):
			if base == search_string[i]:
				score = score
			else:
				score +=1
		score_list.append(score)
# 	else:
# 		score_list = [None, None, None]
	if score_list:
		min_score_index = score_list.index(min(score_list))
	else:
		min_score_index = 0
		score_list = [None]
	
	#returns tuple of search strings starting and ending index
	# and match score
	return (min_score_index, 
			min_score_index + len(search_string), 
			min(score_list))
	
def decide_sequence_direction(sequence, 
							  fwd_tuple1, rev_tuple1, 
							  fwd_tuple2, rev_tuple2):

	#define when to rule out too many errors in binding site
	acceptable_pct_ns_in_binding_site = len(bind1)* 0.25
	
	sequence_direction_list = []
	
	#assigns bind1 direction based on score for scanning
	#binding site 1 forward and reverse
	if (fwd_tuple1[2] < rev_tuple1[2] and 
		fwd_tuple1[2] <= acceptable_pct_ns_in_binding_site 
		and sequence.initial_screen == "OK"):	
		sequence.bind1dir = "forward"
	elif (fwd_tuple1[2] > rev_tuple1[2] and 
		  rev_tuple1[2] <= acceptable_pct_ns_in_binding_site 
		  and sequence.initial_screen == "OK"):
		sequence.bind1dir = "reverse"
	else:
		sequence.bind1dir = "dunno"
		
	
	#assigns bind1 direction based on score for scanning
	#binding site 1 forward and reverse	
	if (fwd_tuple2[2] < rev_tuple2[2] and 
		fwd_tuple2[2] <= acceptable_pct_ns_in_binding_site 
		and sequence.initial_screen == "OK"):	
		sequence.bind2dir = "forward"
	elif (fwd_tuple2[2] > rev_tuple2[2] and 
		  rev_tuple2[2] <= acceptable_pct_ns_in_binding_site 
		  and sequence.initial_screen == "OK"):
		sequence.bind2dir = "reverse"
	else:
		sequence.bind2dir = "dunno"
		
		
	#checks if bind1 and bind2 direction match. 	
	if sequence.bind1dir == sequence.bind2dir:
		sequence.direction = sequence.bind1dir
	else:
		sequence.direction = "conflicting directions"
	
		
		
	

def find_distance_between_bind_sites(sequence, fwd1, rev1, fwd2, rev2):

	if sequence.direction == "forward":
		sequence.distance_between_sites = fwd2[0]-fwd1[1]
		print sequence.name, sequence.distance_between_sites, "forward"
	elif sequence.direction == "reverse":
		sequence.distance_between_sites = rev2[0] - rev1[1]
		print sequence.name, sequence.distance_between_sites, "reverse"
	else:
		print "This sequence didn't have a distance"
		sequence.distance_between_sites = None
	

def main():
	
	#### SETUP ######
	file_to_open = "nrxn1.seq"
	#file_to_open = "epas1b.seq"
	bp_cutoff_value = 275
	n_pct_cutoff_value = 32
	#sets how many digits for numbers defined as Decimal(#)
	getcontext().prec = 3
	
	#running	
	sequence_list = read_fasta(file_to_open)
	
	#test_list = [1,5,10,25,27,28,31,36,39,45,55,56,57,62,64,96]
	#test_list_indexes = [x-1 for x in test_list]
	
	for i, sequence in enumerate(sequence_list):
		#adds basic properties like length, num ns, pct ns
		annotate_sequence(sequence, 
						  bp_cutoff_value, 
						  n_pct_cutoff_value)
		
		#Finds location of best match for binding site 1
		#in forward and reverse directions
		fwd1  = find_best_match_location(sequence, bind1)
		rev1 = find_best_match_location(sequence.reverse_complement(),
									   bind1)
		#Finds location of best match for binding site 2
		#in forward and reverse directions
		fwd2  = find_best_match_location(sequence, bind2)
		rev2 = find_best_match_location(sequence.reverse_complement(),
									   bind2)
		#compares directions for binding site 1 and 2
		#assigns overall sequence direction
		decide_sequence_direction(sequence, 
								  fwd1, rev1, 
								  fwd2, rev2)
		
		find_distance_between_bind_sites(sequence, 
										 fwd1, rev1, 
										 fwd2, rev2)
	
	
main()