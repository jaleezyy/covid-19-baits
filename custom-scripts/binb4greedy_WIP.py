#!/bin/python

### WORK IN PROGRESS
### Goal is to create my own version of the greedy algorithm based on BLAST local alignment
### Run BLASTn on bait set itself
### Isolate a set of bait headers (SeqRecord) - check memory usage
### Parse BLAST file relative to bait header set

import re
import sys, os
from Bio import SeqIO
import csv
from multiprocessing import Pool, Process, Lock, Queue

try:
	input = sys.argv[1]
	blast = sys.argv[2]
	output = sys.argv[3]
	version = int(sys.argv[4])
except IndexError:
	print("""
Arguments are as followed:
	new_greedy.py <input_file> <blast_file> <output_file> <version_1_or_2> <num_threads (default = 1)>
	
	For <version> type either "1" or "2" or "3":
		Version "1" Greedy Algorithm first analyses all alignment hits for a bait prior to removal.
		Version "2" Greedier Algorithm removes baits by analysing each BLAST line once (line-by-line removal).
		""")
	exit()

bait_set = {}
for seq_record in SeqIO.parse(input, 'fasta'): # extract all baits - subject to removal)
	bait_set[seq_record.description] = seq_record.seq
	
print("Bait Set Before: \n")
print(bait_set)
print("\n")


if version == 1:	
	print("Running greedy algorithm!\n")
	analysed = set()
	with open(blast) as b:
		for line in b:
			if str(line).split(",")[0] not in bait_set: # if already removed, no need to analyse it
				continue
			else:
				hits = []
				bait = str(line.split(",")[0])
				if bait in analysed:
					continue
				else:
				#print(bait)
					for l in open(blast):
						if l.startswith(bait) is True:
							hits.append(l)
						else:
							pass
					analysed.add(bait) # want to avoid duplicate analyses
								
		#print(hits)
			if len(hits) >= 1:
				for i in hits:
					if bait.split("_")[0] in str(i).split(",")[1]:
						continue
					else:
						length = float(i.split(",")[3])
						similarity = float(i.split(",")[2])
						if similarity >= float(80):
							if length >= float(60):
								try:
									del bait_set[i.split(",")[1]]
								except KeyError:
									print("Bait " + str(i.split(",")[1]) + " not found, skipping.")
									continue
							else:
								continue
						else:
							continue
			else:
				continue

elif version == 2:
	print("Running greedier algorithm!\n")
	with open(blast) as b:
		for line in b:
			bait = str(line).split(",")[0]
			if bait not in bait_set: # if already removed
				continue
			elif bait.split("_")[0] in str(line).split(",")[1]:
				continue
			else:
				length = float(line.split(",")[3])
				similarity = float(line.split(",")[2])
				if similarity >= float(80):
					if length >= float(60):
						try:
							del bait_set[line.split(",")[1]]
						except KeyError:
							print("Bait " + str(line.split(",")[1]) + " not found, skipping.")
							continue
					else:
						continue
				else:
					continue
else:
	exit("Incorrect version number in position 4!")

print("\nBait Set After: \n")
print(bait_set)
print("\n")

with open(output, 'w+') as out:
	for bait in bait_set:
		out.write(">" + str(bait) + "\n" + str(bait_set[bait]) + "\n")
		
print("Done.")
		
			
			
		
		
	
	
	
