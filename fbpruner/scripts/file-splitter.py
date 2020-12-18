 #!/usr/bin/env python
###############################################################################
#  	   FILE: file_splitter.py
# 	PURPOSE: Takes a file and splits it by number of lines.
#
#  	 AUTHOR: Dave Rich
#  	    BUG: Lots.
###############################################################################

import sys

def split_fasta( fname, num_per_split ):
	# create template for outnames
	# special case if filename has no extension
	if fname.find(".") == -1:
		file = fname + "."
		ext = ""
	# otherwise
	else:
		fnames = fname.split(".")
		file = ""
		ext = "." + fnames[len(fnames)-1]
		for i in range(len(fnames)-1):
			file += fnames[i] + "."

	# line count
	l_cnt = 0
	# model count
	m_cnt = 0
	# file count
	f_cnt = 0
	# outfile name
	outname = f'{file}{f_cnt}{ext}'
	fout = open(outname, "w+")
	print(f"Adding file: {outname}...")

	with open(fname, "r") as fp:
		for line in fp:
			# indicates the start of a new fasta entry
			if line.startswith(">"):
				m_cnt += 1
				# if we've reached the number per file, create new file
				if m_cnt >= num_per_split:
					m_cnt = 0
					f_cnt += 1
					# swap outfiles
					fout.close()
					outname = f"{file}{f_cnt}{ext}"
					fout = open(outname, "w+")
					print(f"Adding file: {outname}...")

			# add line to current file
			fout.write(line)

	# close final file
	fout.close()

	
def split_hmm( fname, num_per_split ):
	# create template for outnames
	# special case if filename has no extension
	if fname.find(".") == -1:
		file = fname + "."
		ext = ""
	# otherwise
	else:
		fnames = fname.split(".")
		file = ""
		ext = "." + fnames[len(fnames)-1]
		for i in range(len(fnames)-1):
			file += fnames[i] + "."

	# line count
	l_cnt = 0
	# model count
	m_cnt = 0
	# file count
	f_cnt = 0
	# outfile name
	outname = f'{file}{f_cnt}{ext}'
	fout = open(outname, "w+")
	print(f"Adding file: {outname}...")

	with open(fname, "r") as fp:
		for line in fp:
			# add line to current file
			fout.write(line)

			# indicates the start of a new fasta entry
			if line.startswith("//"):
				m_cnt += 1
				# if we've reached the number per file, create new file
				if m_cnt >= num_per_split:
					m_cnt = 0
					f_cnt += 1
					# swap outfiles
					fout.close()
					outname = f"{file}{f_cnt}{ext}"
					fout = open(outname, "w+")
					print(f"Adding file: {outname}...")

			
	# close final file
	fout.close()


##############################################################################
###########################         MAIN         #############################
##############################################################################

if len(sys.argv) != 4:
	print("Usage: [hmm/fasta] <file> <num_of_models_per_file>")
	exit(0)

if sys.argv[1] == "fasta":
	ftype = "fasta"
elif sys.argv[1] == "hmm":
	ftype = "hmm"
else:
	print(f'"{sys.argv[1]} is not a valid file type. Must be one of [hmm/fasta]."')
	exit(0)

fname = sys.argv[2]
num_per_split = int(sys.argv[3])

print(f"FILE: {fname}, FILETYPE: {ftype}, NUM_PER_SPLIT: {num_per_split}")

split_fasta( fname, num_per_split )
