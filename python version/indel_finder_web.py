from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from decimal import *
import csv, re

def analyze_wt(wt, b1, b2):
	b1_end_ind = wt.seq.find(b1) + len(b1)
	b2_start_ind = wt.seq.find(b2)
	wt_distance_between_sites = b2_start_ind - b1_end_ind
	return wt_distance_between_sites

def create_temp_file(sequences):
	temp_file = open("temp.seq", 'w+b')
	temp_file.write(sequences)
	temp_file.close()

	return "temp.seq"

def read_fasta(file_to_open):	
	seq_file = open(file_to_open, "rU")
	sequence_list = []
	
	for record in SeqIO.parse(seq_file, "fasta", 
							  alphabet = generic_dna):
		sequence_list.append(record)
		record.initial_screen = "Omitted" 
		record.length = None 
		record.n_count = None 
		record.n_pct = None 
		record.direction = None 
		record.bind1dir = None 
		record.bind2dir = None 
		record.distance_between_sites = None
		record.score = None
		record.b1 = None
		record.b2 = None
		
	seq_file.close()
	return sequence_list
	
def annotate_sequence(sequence, bpcutoff, npctcutoff):
	#sets how many digits for numbers defined as Decimal(#)
	getcontext().prec = 3
	
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

	if score_list:
		min_score_index = score_list.index(min(score_list))
	else:
		min_score_index = 0
		score_list = [None]
	
	#returns tuple of search strings starting and ending index
	#and match score
	return (min_score_index, 
			min_score_index + len(search_string), 
			min(score_list))

	
def decide_sequence_direction_for_one_only(sequence,
 										   fwd_tuple,
 										   rev_tuple,
 										   bind1):

	#define when to rule out too many errors in binding site
	acceptable_pct_ns_in_binding_site = len(bind1)* 0.25

 	if (fwd_tuple[2] < rev_tuple[2] and 
		fwd_tuple[2] <= acceptable_pct_ns_in_binding_site 
		and sequence.initial_screen == "OK"):	
		return "forward"
	elif (fwd_tuple[2] > rev_tuple[2] and 
		  rev_tuple[2] <= acceptable_pct_ns_in_binding_site 
		  and sequence.initial_screen == "OK"):
		return "reverse"
	else:
 		return "dunno"  

										
def decide_sequence_direction(sequence, 
							  fwd_tuple1, rev_tuple1, 
							  fwd_tuple2, rev_tuple2,
							  bind1):
	
 	sequence.bind1dir = decide_sequence_direction_for_one_only(
 		sequence, fwd_tuple1, rev_tuple1, bind1)
 	sequence.bind2dir = decide_sequence_direction_for_one_only(
 		sequence, fwd_tuple2, rev_tuple2, bind1)
		
	#checks if bind1 and bind2 direction match. 	
	if sequence.bind1dir == sequence.bind2dir:
		sequence.direction = sequence.bind1dir
	else:
		sequence.direction = "conflicting directions"
	

def find_distance_between_bind_sites(sequence, fwd1, rev1, fwd2, rev2):

	if sequence.direction == "forward":
		sequence.distance_between_sites = fwd2[0]-fwd1[1]
		sequence.score = (fwd1[2], fwd2[2])
		sequence.b1 = sequence.seq[fwd1[0]:fwd1[1]]
		sequence.b2 = sequence.seq[fwd2[0]:fwd2[1]]
	elif sequence.direction == "reverse":
		sequence.distance_between_sites = rev2[0] - rev1[1]
		sequence.score = (rev1[2], rev2[2])
		sequence.b1 = sequence.seq.reverse_complement()[rev1[0]:rev1[1]]
		sequence.b2 = sequence.seq.reverse_complement()[rev2[0]:rev2[1]]
	else:
		sequence.distance_between_sites = None
	

def make_html_table(sequence_list):
	#rdr = csv.reader(open("summary.csv", "r"))
	csv_data = []
	csv_data.append(["well name", "SPACER LENGTH", "initial screen", "seq length", "pct n's",
					 "seq direction", "bind seq 1", "bind seq 2", "match score"])
	for well in sequence_list:

		if well.direction == "forward" or well.direction == "reverse":
			csv_data.append([well.name, well.distance_between_sites, well.initial_screen, well.length,
			 			 	 well.n_pct, well.direction, well.b1, well.b2, well.score])
		else:
			csv_data.append([well.name, well.distance_between_sites, well.initial_screen, well.length,
						 	 well.n_pct, well.direction, "none", "none", "none"])

	return csv_data

def write_fasta(output_list, wt):
	seqs_same_dir = [wt]
	for seq in output_list:
		if seq.direction == "reverse":
			seq.seq = seq.seq.reverse_complement()
			seqs_same_dir.append(seq)			
		else:
			seqs_same_dir.append(seq)
	
	output_handle = open("align.seq", "w+b")
	SeqIO.write(seqs_same_dir, output_handle, "fasta")
	output_handle.close()


def main_for_loop(sequence_list, bind1, bind2, wt_dist):

	bp_cutoff_value = 275
	n_pct_cutoff_value = 32
	output_list = []

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
								  fwd2, rev2,
								  bind1)
		
		find_distance_between_bind_sites(sequence, 
										 fwd1, rev1, 
										 fwd2, rev2)

		#makes list to put into fasta file
		if (sequence.distance_between_sites != wt_dist and
			sequence.distance_between_sites != None):
			output_list.append(sequence)

	return output_list


def sanitize_wt(wt):
	wt = wt.replace(" ", "")
	return wt

def make_string_seqrecord(wt, id_txt, name = " ", desc = " " ):

	wt = SeqRecord(Seq(wt, generic_dna), id = "wt", name = "wt", description = "wt")
	#wt.id = wt
	wt.direction = "forward"
	return wt

def main(b1, b2, gene_name, wt, sequences):
	
	#### SETUP ######
	bind1 = b1
	bind2 = b2
	sequences = sequences
	gene_name = gene_name
	wt = sanitize_wt(wt)
	wt = make_string_seqrecord(wt, "wt")
	wt_dist = analyze_wt(wt, b1, b2)
	file_to_open = create_temp_file(sequences)
	sequence_list = read_fasta(file_to_open)
	output_list = main_for_loop(sequence_list, bind1, bind2, wt_dist)
	html_data = make_html_table(sequence_list)
	write_fasta(output_list, wt)
	return html_data
