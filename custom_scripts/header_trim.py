#!/bin/python

### Goal is to trim header and output FASTA file with headers modified
### Header is modified to be AccessionID_<bait_range> i.e. A10048.01_21-100

from Bio import SeqIO
import os, sys
from multiprocessing import Pool, Process, Lock, Queue

try:
	input = sys.argv[1]
	threads = int(sys.argv[2])
except IndexError:
	try:
		input = sys.argv[1]
		threads = 1
	except IndexError:
		print("""
Arguments are as followed:
	header_trim.py <input_file> <num_threads (optional)>
	""")
		exit()


### https://biopython.org/wiki/Split_large_file - iterator for FASTA file

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
	num = len([1 for line in open(input) if line.startswith(">")])
	return num

def split_fasta(fasta_file, num, threads):
	batch_all = [] # store split input
	record_iter = SeqIO.parse(fasta_file,"fasta")
	if threads <= 0:
		raise ValueError
	for i, batch in enumerate(batch_iterator(record_iter, int(num//threads + (num % threads > 0)))):
		i = i + 1
		if i <= threads:
			batch_all.append(batch)
		else:
			break
	return batch_all


def trim_header(batch):
	header = []
	seq = []
	for seq_record in batch:
		header.append(str(str(seq_record.description).split(" ")[0] + "_" + str(seq_record.description).rsplit("_",1)[1]))
		seq.append(str(seq_record.seq))
	
	return header,seq
	
def run():
	num_seq = count_seq(input)
	batch = split_fasta(input, num_seq, threads)
	print("Split into " + str(len(batch)) + " subunit(s).")
	p = Pool(threads)
	with open(input + "_trimmed", "w") as out:
		for i,j in p.imap_unordered(trim_header, batch):
			for h,s in zip(i,j):
				out.write(">" + str(h) + "\n" + str(s) + "\n")

if __name__ == '__main__':
	run()
	print("Trimmed headers.")
	exit()
