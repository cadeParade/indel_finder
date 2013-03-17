from Bio import SeqIO
from Bio.Alphabet import generic_dna
import numpy as np
import string
import re
from decimal import *

def read_fasta(file_to_open):
	
	seq_file = open(file_to_open, "rU")
	sequence_list = []
	
	for record in SeqIO.parse(seq_file, "fasta", 
							  alphabet = generic_dna):
		sequence_list.append(record)
		
	seq_file.close()
	return sequence_list
	
def annotate_sequences(sequence_list, bpcutoff, npctcutoff):
	
	sequence_dict = {}
	
	for i,sequence in enumerate(sequence_list):
	
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
	
	
	print sequence_list[1].name, "name"
	print sequence_list[1].length, "length"
	print sequence_list[1].n_count, "n count"
	print sequence_list[1].n_pct, "n pct"
	print sequence_list[1].initial_screen, "initial screen"
	print sequence_list[1]
	

def main():
	
	#### SETUP ######
	file_to_open = "epas1b.seq"
	bp_cutoff_value = 275
	n_pct_cutoff_value = 32
	#sets how many digits for numbers defined as Decimal(#)
	getcontext().prec = 3
	
	
	#running	
	sequence_list = read_fasta(file_to_open)
	annotate_sequences(sequence_list,
					   bp_cutoff_value, 
					   n_pct_cutoff_value)
	
main()