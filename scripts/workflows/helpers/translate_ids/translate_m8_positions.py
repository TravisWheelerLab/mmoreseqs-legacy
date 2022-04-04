#!/bin/usr/env python
#####################################################################
#  - FILE:  translate_m8_positions.py
#  - DESC:   Translate .m8 alignments from positions in MMseqs models ("p") to the corresponding
# 			 positions in the HMMER models, based on HMMER consensus sequence ("s").
# 			 
#####################################################################

import os,sys
import numpy as np

debug = True

# debug print
def dprint(*args):
	if (debug):
		print(args)
	return 

# load file with identity alignment b/t mmseqs and hmm models.
def load_identity(my_path, my_data):
	my_fp = open(my_path, "r")

	line = my_fp.readline()
	while line:
		fields 	= line.split()
		p_name	= fields[0]
		# cigar 	= fields[len(fields)-1]

		# parse data entry
		my_data[p_name] = {}
		my_data[p_name]["general"] = fields
		p_name 							= fields[0]
		s_name 							= fields[1]
		p_len 							= int(fields[2])
		p_beg 							= int(fields[3])
		p_end 							= int(fields[4])
		p_len 							= int(fields[2])
		s_beg 							= int(fields[6])
		s_end 							= int(fields[7])
		alnlen 							= int(fields[8])
		cigar 							= fields[9]

		# convert cigar alignment to array of p_positions and s_positions
		p_pos, s_pos = cigar_to_array(cigar, p_beg, p_end, s_beg, s_end)
		# my_data[p_name]["states"]  = states 
		# my_data[p_name]["runs"]  	= runs 
		my_data[p_name]["p_pos"] 	= p_pos
		my_data[p_name]["s_pos"] 	= s_pos

		line = my_fp.readline()

	my_fp.close()
	return 

# convert cigar alignment to array of index p_positions and s_positions
def cigar_to_array(cigar, p_beg, p_end, s_beg, s_end):
	# append to begining of list positions before the alignment
	p_pos = [0, p_beg]
	s_pos = [0, s_beg]

	# break up cigar alignment... 
	# ...into lists of states
	states = cigar
	# remove all digits
	for i in range(10):
		states = states.replace(str(i), " ")
	states = states.split()
	# ...into list of runs
	runs = cigar
	# remove all state letters (M,I,D)
	runs = runs.replace("M"," ").replace("I"," ").replace("D"," ")
	runs = runs.split()
	runs = [int(x) for x in runs]
	# ...into list of positions
	p_cnt = p_beg
	s_cnt = s_beg
	for i in range(len(runs)):
		if (states[i] == "M"):
			p_cnt += runs[i]
			s_cnt += runs[i]
		if (states[i] == "I"):
			p_cnt += runs[i]
			s_cnt += 0
		if (states[i] == "D"):
			p_cnt += 0
			s_cnt += runs[i]
		p_pos.append(p_cnt)
		s_pos.append(s_cnt)
	
	# append to end of list positions after the alignment 
	p_pos.append(p_end)
	s_pos.append(s_end)
	
	# print("P:({},{}), S:({},{})".format(p_beg,p_end,s_beg,s_end))
	# print("CIGAR:", cigar)
	# print("RUNS:", len(runs), runs)
	# print("STATES:", len(states), states)

	# print("P_POS:", len(p_pos), p_pos)
	# print("S_POS:", len(s_pos), s_pos)

	return p_pos, s_pos

# translate position of points b/t mmseqs model position -> hmmer model position
def translate_seed_range(fields, id_fields, my_stats):
	# parse data from mmseqs alignments
	p_name 		= fields[0]
	p_alnbeg 	= int(fields[6])
	p_alnend 	= int(fields[7])
	# print("P_NAME:", p_name)
	# print("ID_FIELDS:", id_fields)

	# parse data from model identity alignments
	p_len 		= int(id_fields["general"][2])
	p_beg 		= int(id_fields["general"][3])
	p_end			= int(id_fields["general"][4])
	s_len 		= int(id_fields["general"][5])
	s_beg 		= int(id_fields["general"][6])
	s_end 		= int(id_fields["general"][7])
	p_pos 		= id_fields["p_pos"]
	s_pos 		= id_fields["s_pos"]

	# output
	s_alnbeg 	= -1
	s_alnend 	= -1

	# translate each position for p_pos -> s_pos
	my_translate = translate_seed_position_1
	# print("P_POS:", len(p_pos), p_pos)
	# print("S_POS:", len(s_pos), s_pos)
	s_alnbeg, stat_beg = my_translate(p_alnbeg, s_pos, s_beg, s_end, s_len, p_pos, p_beg, p_end, p_len, my_stats)
	s_alnend, stat_end = my_translate(p_alnend, s_pos, s_beg, s_end, s_len, p_pos, p_beg, p_end, p_len, my_stats)
	my_stats[(stat_beg,stat_end)] += 1
	# print("P:({},{}) -> S:({},{})".format(p_alnbeg, p_alnend, s_alnbeg, s_alnend))
	# print("STATS:", stat_beg, stat_end)

	return s_alnbeg, s_alnend

# Method #1: translate position of points b/t mmseqs model position -> hmmer model position
#  - 	This method handles positions outside of the alignment by approximating the location in the
# 		hmmer model by the percentage of the outer nodes outside the alignment in the mmseqs model.
def translate_seed_position_1(p_alnpt, s_pos, s_beg, s_end, s_len, p_pos, p_beg, p_end, p_len, my_stats):
	# output
	s_alnpt = -1

	# if alignment is before identity alignment 
	if ( p_alnpt < p_beg ):
		my_stat = "before_aln"
		p_perc = 0
		if ( p_beg > 1 ):
			p_perc 	= float(p_beg - p_alnpt) / float(p_beg)
		s_alnpt = int(s_beg - ( p_perc * float(s_beg) ) )
		s_alnpt = max(0, s_alnpt)
		pass

	# if alignment is inside identity alignment 
	elif ( p_alnpt <= p_end ):
		my_stat = "inside_aln" 
		for i in range(1,len(p_pos)):
			if ( p_pos[i] >= p_alnpt ):
				if ( s_pos[i] != s_pos[i-1] ):
					s_alnpt = s_pos[i-1] + (p_alnpt - p_pos[i-1])
				else:
					s_alnpt = s_pos[i-1]
				break
		pass

	# if alignment is after identity alignment 
	elif ( p_alnpt > p_end ):
		my_stat = "after_aln"
		p_perc = 0
		if ( p_len - p_end != 0 ):
			p_perc 	= float(p_alnpt - p_end) / float(p_len - p_end)
		s_alnpt = int( s_end + ( p_perc * float(s_len - s_end) ) )
		s_alnpt = min(s_len-1, s_alnpt)
		pass 

	# shouldn't get here
	else:
		pass
	return s_alnpt, my_stat

# Method #2: translate position of points b/t mmseqs model position -> hmmer model position
# 	- 	This method handles positions outside of the alignment by truncating the results  
# 		positions within the model alignment, then truncating the target sequence accordingly
# 		(required that mmseqs search for query/target include full traceback alignment "-a")
def translate_seed_position_2(p_alnpt, s_pos, s_beg, s_end, s_len, p_pos, p_beg, p_end, p_len, my_stats):
	s_alnpt = -1
	# if alignment is before identity alignment 
	if ( p_alnpt < p_beg ):
		my_stat = "before_aln"
		p_perc = 0
		if ( p_beg > 1 ):
			p_perc 	= float(p_beg - p_alnpt) / float(p_beg - 1)
		s_alnpt = s_beg - int( p_perc * float(s_beg) )
		pass
	# if alignment is inside identity alignment 
	elif ( p_alnpt <= p_end ):
		my_stat = "inside_aln" 
		for i in range(1,len(p_pos)):
			if ( p_pos[i] >= p_alnpt ):
				if ( s_pos[i] != s_pos[i-1] ):
					s_alnpt = s_pos[i-1] + (p_alnpt - p_pos[i-1])
				else:
					s_alnpt = s_pos[i-1]
				break
		pass
	# if alignment is after identity alignment 
	elif ( p_alnpt > p_beg ):
		my_stat = "after_aln"
		p_perc = 0
		if ( p_len - p_end != 0 ):
			p_perc 	= float(p_alnpt - p_end) / float(p_len - p_end)
		s_alnpt = s_end + int( p_perc * float(s_len - s_end) )
		pass 
	# shouldn't get here
	else:
		pass
	return s_alnpt, my_stat

#####################################################################
##  MAIN  ###########################################################
#####################################################################

my_args 	= {}
my_fps	= {}
my_data 	= {}
my_stats = {}

my_data["res_m8"] = {}
my_data["m2m_m8"] = {}
my_data["out_m8"] = {}

req_args = 3
if len(sys.argv) < req_args + 1:
	print("# ERROR: incorrect number of args.")
	print("# ./ <i:mmseqs_model_results_m8> <i:model2model_m8> <o:hmmer_model_results_m8>")
	exit(-1)

# load args
my_args["res_m8"] = sys.argv[1]
my_args["m2m_m8"] = sys.argv[2]
my_args["out_m8"] = sys.argv[3]

# load identity data
load_identity(my_args["m2m_m8"], my_data["m2m_m8"])

# translate alignment seeds
my_fps["res_m8"] = open(my_args["res_m8"], "r")
my_fps["out_m8"] = open(my_args["out_m8"], "w")

for pos_beg in ["before_aln", "inside_aln", "after_aln"]:
	for pos_end in ["before_aln", "inside_aln", "after_aln"]:
		my_stats[(pos_beg,pos_end)] = 0

line = my_fps["res_m8"].readline()
while line:
	try:
		# get input data from mmseqs result entry
		fields 		= line.split()
		p_name 		= fields[0]
		# get input data from model identity alignments
		id_fields 	= my_data["m2m_m8"][p_name]
		# translate viterbi seed positions to hmm model positions
		s_alnbeg, s_alnend = translate_seed_range(fields, id_fields, my_stats)

		# replace mmseqs profile indexes with hmmer profile indexes
		# print("FIELDS:", len(fields), fields)
		entry = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format( 
			fields[0], fields[1], fields[2], fields[3], fields[4],  fields[5], 
			int(s_alnbeg),  int(s_alnend),  fields[8], fields[9], fields[10], fields[11] )
		my_fps["out_m8"].write( "{}\n".format(entry) )

	# catch keyerror
	except:
		e = sys.exc_info()[0]
		print("ERROR: %s" % (e) )
	line = my_fps["res_m8"].readline()

my_fps["res_m8"].close()
my_fps["out_m8"].close()

print("# completed successfully.")

pos = ["before_aln", "inside_aln", "after_aln"]
for pos_beg in pos:
	for pos_end in pos:
		print("# stats({},{}): {}".format(pos_beg, pos_end, my_stats[(pos_beg,pos_end)]) )
