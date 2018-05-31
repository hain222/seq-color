#!/usr/bin/env python3

# -------------------------------------------------------------------

# File: seqColor.py
# Author: Harrison Inocencio
# Date: 05-29-18
# Purpose: Color codes sequences based on their quality scores. 

# Usage: NA

# Compilation: NA

# Notes:
# 1. Should use some kind of indexing algorithm (fast) to search for a fasta entries sister fastq
#	 entry. Needs more thought.
# 2. 
# 3.
# 4.
# 5.

# TODO:
# 1. Update Help messages
# 2. Needs Sequence/id search function
# 3.
# 4.
# 5.

# -------------------------------------------------------------------

import os
import argparse
import colorama
from Bio import SeqIO

# Globals:
wrap_count = 70

bad_erange = 19

avg_srange = 20
avg_erange = 28

good_srange = 29

# crit_error function
# > Handles unrecoverable program errors
def crit_error(error):
	print("ERROR:", error)
	exit(1)

# fasta_print function
# > If -f is specified, prints all records to stdout in fasta format
def fasta_print(rec_list):
	global wrap_count	

	try:
		for rec in rec_list:
			wrap_iter = 0
			quals = rec.letter_annotations["phred_quality"]
			print(">"+rec.id)
		
			for idx, base in enumerate(rec.seq):
				wrap_iter += 1
				qual_base_print(quals[idx], base)
				if wrap_iter == wrap_count:
					print()
					wrap_iter = 0
			print()
	except KeyboardInterrupt:
		print()
		exit(0)

	return 0

# seq_iterator function
# > iterates over the different entries, printing them one at a time
def seq_iterator(rec_list):
	
	try:
		for rec in rec_list:	
			seq_print(rec)
			input("\nPress a key to continue...")
		print("All Sequences Processed...")
	except KeyboardInterrupt:
		print()
		exit(0)

# qual_base_print function
# > prints a base according to an associated quality score and returns
# > the given score numerically
def qual_base_print(qual_score, base):
	global bad_erange
	global avg_srange
	global avg_erange
	global good_srange

	colorama.init()

	# Calls
	# 0 = initial
	# 1 = bad
	# 2 = average
	# 3 = good
	call = 0;

	if qual_score <= bad_erange:
		print(colorama.Back.RED + base, end='')
		print(colorama.Style.RESET_ALL, end='')
		call = 1
	elif avg_srange <= qual_score <= avg_erange:
		print(colorama.Back.YELLOW + base, end='')
		print(colorama.Style.RESET_ALL, end='')
		call = 2
	elif qual_score >= good_srange:
		print(colorama.Back.GREEN + base, end='')
		print(colorama.Style.RESET_ALL, end='')
		call = 3

	return call


# seq_print function
# > reprints the given seq record according to the colorcode
def seq_print(rec):
	global wrap_count

	wrap_iter = 0
	bad_count = 0
	quals = rec.letter_annotations["phred_quality"]
	
	print("######################################################################\n\n"+rec.id)
	for idx, base in enumerate(rec.seq):
		wrap_iter += 1

		call = qual_base_print(quals[idx], base)
		if call == 1:
			bad_count += 1

		if wrap_iter == wrap_count:
			print()
			wrap_iter = 0

	per_bad = round((bad_count/len(rec.seq) * 100), 1)
	print("\n\n% Bad: "+str(per_bad))
	print("Read Length: "+str(len(rec.seq)))
	print("######################################################################")
			

# main function
# > Parses input and executes primary program functions
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('fastq_input', help='fastq_input HELP')
	parser.add_argument("-f", "--fasta", help="--fasta HELP", action="store_true")
	parser.add_argument("-m", "--based-on-motif", help="--based-on-motif HELP", action="store_true")
	args = parser.parse_args()

	fastq_fid = args.fastq_input
	seq_recs = []

	try:
		for seq_record in SeqIO.parse(fastq_fid, "fastq"):
			seq_recs.append(seq_record)
	except FileNotFoundError as err:
		crit_error(err)


	#print(seq_recs[0].letter_annotations["phred_quality"])
	
	#print(args.based-on-motif)

	if args.fasta:
		fasta_print(seq_recs)
	else:
		seq_iterator(seq_recs)

	return 0
main()
