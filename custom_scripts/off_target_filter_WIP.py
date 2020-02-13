#!/bin/python

### WORK IN PROGRESS

from Bio import SeqIO
import argparse
import re
import os, sys
import glob
import csv
import datetime
import textwrap
from multiprocess import Pool, Process, Lock, Queue

### Goal is to filter bait sequences based on BLASTn alignment (with multiprocessing)
### Want to remove BLASTn results linked to Viruses first
### Remaining results are to be used to remove baits with matching elements in the header (i.e. AccID_rng)

### BLASTn should be run with the following: -outfmt '10 std staxids sscinames sblastnames sskingdoms'

### Not meant to function as filter within the keep_term (i.e. will not distinguish hybridization within Viruses)

def create_parser():
	parser = argparse.ArgumentParser(description='Filter bait sequences based on BLAST alignment.')
	parser.add_argument('-i', '--input', dest="input_fasta_file", required=True, help="input FASTA file")
	parser.add_argument('-b', '--blast', dest="input_blast_file", required=True, help="input BLAST results file")
	parser.add_argument('-o', '--output', dest="output_fasta_file", required=True, help="output fasta file")
	parser.add_argument('-t', '--threads', dest="num_threads", required=False, help="number of threads - default is 1")
	parser.add_argument('-k', '--keep', dest="keep_term", required=False, help="keep baits with matching sskingdoms only - default filters based on all information without exception.")
	return parser

def batch_iterator(iterator, batch_size): # create iterator
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = next(iterator)
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch

def count_seq(input):
	num_seq = len([1 for line in open(input) if line.startswith(">")])
	return num_seq
	
def count_line(input):
	num_line = len([1 for line in open(input)])
	return num_line

def generator_file(file):
	with open(file) as r:
		for line in r:
			yield line

def split_file(file, num, threads):
	batch_all = [] # store split input
	if ".fasta" in file:
		record_iter = SeqIO.parse(file,"fasta")
	elif ".fa" in file:
		record_iter = SeqIO.parse(file,"fasta")
	elif ".fa_trimmed" in file:
		record_iter = SeqIO.parse(file,"fasta")
	elif ".fasta_trimmed" in file:
		record_iter = SeqIO.parse(file,"fasta")
	else:
		record_iter = iter(generator_file(file))
	if threads <= 0:
		raise ValueError
	for i, batch in enumerate(batch_iterator(record_iter, int(num//threads + (num % threads > 0)))):
		i = i + 1
		if i <= threads:
			batch_all.append(batch)
		else:
			break
	return batch_all
		

def parse_blast(blast, keep_term): # function acts to remove entries containing the sskingdom term (will have to expand)
	keep_list = [] # keep track of sequences related with keep_term
	blast_list = []
	for batch in blast:
		for line in batch:
			if str(keep_term).upper() in str(line).upper().rsplit(",",1)[1]:
				if str(line).split(",")[0] in keep_list:
					pass
				else:
					keep_list.append(str(line).split(",")[0]) # used to check bait header
			else:
				blast_list.append(line) #full blast line --> to tmp file
	return keep_list, blast_list

def to_remove(seq, blast): # generates list of sequence headers that are NOT to be found in the final FASTA file
	input=seq
	removal_list = []
	with open(blast) as b:
		for line in b:
			ref = line.split("\n")[0]
			#print(ref)
			for ele in input:
				for seq_record in ele:
					query_1 = str(seq_record.description).split(" ")[0]
					query_2 = str(seq_record.description).rsplit("_",1)[1]
					if query_1 in str(ref):
						if query_2 in str(ref):
							removal_list.append(str(query_1+"_"+query_2))
							break
						else:
							pass
					else:
						pass
	return removal_list			

def off_target(seq, removal_list, blast): # filters the input FASTA, outputting new FASTA file reduced based on hits found in BLAST output
	input=seq
	off_list = []
	for ele in seq:
		for seq_record in ele:
			query_1 = str(seq_record.description).split(" ")[0]
			query_2 = str(seq_record.description).rsplit("_")[1]
			if str(query_1+"_"+query_2) in removal_list:
				with open(blast, 'r') as b:
					for line in b:
						if str(query_1+"_"+query_2) in line.split("\n")[0]:
							if float(str(line).split(",")[3]) > 50.0:
								if float(str(line).split(",")[2]) > 80.0:
									continue
								else:
									off_list.append((str(seq_record.description), str(seq_record.seq)))
									#out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")
							else:
								off_list.append((str(seq_record.description), str(seq_record.seq)))
								#out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")
						else:
							pass
				#continue
			else:
				off_list.append((str(seq_record.description), str(seq_record.seq)))
				#out.write(">" + str(seq_record.description) + "\n" + str(seq_record.seq) + "\n")
	return off_list			
			
					
		#print(str(seq_record.description).split(" ")[0]+"_"+str(seq_record.description).rsplit("_",1)[1])
		#exit()
	
	
def run():
	parser = create_parser()
	args = parser.parse_args()
	if args.num_threads is not None:
		threads = int(args.num_threads)
	else:
		threads = 1
		
	print("Parsing BLAST output.\n")
	
	num_blast = count_line(args.input_blast_file)
	batch_blast = split_file(args.input_blast_file, num_blast, threads)
	print("\nSplit BLAST input into " + str(len(batch_blast)) + " subunit(s).\n")
	#print(batch_blast)
	
	num_seq = count_seq(args.input_fasta_file)
	batch_seq = split_file(args.input_fasta_file, num_seq, threads)
	print("\nSplit sequence input into " + str(len(batch_seq)) + " subunit(s).\n")
	#print(batch_seq)

### Parse BLAST output file

	if args.keep_term is not None:
		print("\nKeeping sequences belonging to sskingdom: " + str(args.keep_term) + "\n")
		keepage = []
		p = Pool(threads)
		with open(args.input_blast_file + "_tmp", "w") as out:
			for i,j in p.starmap(parse_blast, [(batch_blast, args.keep_term)]):
				for line in j: ### fix nested for loop
					out.write(line)
				keepage = i
		blast=args.input_blast_file+"_tmp"
	else:
		p = Pool(threads)
		print("\nNo parsing of BLAST output.\n")
		keepage = []
		blast=args.input_blast_file
	
### Using BLAST output file (parsed or otherwise), figure out which baits need to be removed	

	print("\nDetermining baits to remove.\n")
	
	for i in p.starmap(to_remove, [(batch_seq, blast)]):
		removal = i

	#print("Removed: " + str(removal))
	#print("Kept: " + str(keepage))
	
### Off-target filter using above list of baits to remove (or ignore when checking each bait)	
	
	print("\nFiltering baits.\n")
	count_off_target = 0
	with open(args.output_fasta_file, "w") as out2:
		for i in p.starmap(off_target, [(batch_seq, removal, blast)]):
			for t in i:
				out2.write(">" + str(t[0]) + "\n" + str(t[1]) + "\n")
				count_off_target = count_off_target + 1
	
	count_off_target = num_seq - count_off_target

	return (count_off_target)

if __name__ == '__main__':
	start=str(datetime.datetime.utcnow()).split('.')[0] # start of program
	num_baits = run()
	finish=str(datetime.datetime.utcnow()).split('.')[0] # end of program
	print("\nRemoved " + str(num_baits) + " bait(s).\n")
	print("Start: " + start)
	print("Finish: " + finish)
	exit()