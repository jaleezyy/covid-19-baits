#!/bin/python

### WORK IN PROGRESS

import os, argparse 
import subprocess
import sys
import re
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

# melt_parse.sh executes melt.pl with subsequence parsing (MUST BE IN SAME PATH AS THIS PYTHON SCRIPT)
# After running melt.pl determine which probes are bad relative to melting temperature (Tm)
# Takes 2 files: probes and melting point (mp) information
# Will compare each value in mp and if greater than specified threshold then it will take the header + sequence from probes to new filtered file

userParameters = sys.argv[1] # input FASTA file
low_temp = sys.argv[2] # input minimum temperature
high_temp = sys.argv[3] # input maximum temperature

# if statement determining whether Tm data already exists as proper name --> otherwise run melt_parse.sh
# allow for user to skip melt_parse.sh if already done

if os.path.exists("./Tm_%s.txt" %(userParameters)):
    print("Melting temperatre data exists for original FASTA file, skipping melt_parsh.sh" + "\n")
    Tm_userParameters = "Tm_"+userParameters+".txt"
else:
    print("Running melt.pl on given probe set")
    subprocess.run("~/jalees/Scripts/melt_parse.sh -s %s" %(userParameters), check=True, shell=True)
    print("Generated list of melting temperatures" + "\n") # Tm_userParameters.txt
    Tm_userParameters = "Tm_"+userParameters+".txt" 


print("Original dataset:")
print(userParameters + "\n") # original fasta file with probes


with open(Tm_userParameters, "r") as f:
    lines = f.read().splitlines()
    

print("Melting temperature data:")
print(Tm_userParameters + "\n")
#print(len(lines)) # check length of the list from Tm data which should match length of original fasta file
#print(lines[0]) # check if first value shows

def melt_filter_range(fasta_file, low_temp, high_temp):
    # Create our hash table to add the sequences
    sequences={} # set object to eliminate duplicates
    no_seq={} 
    i = 0

    if isinstance(low_temp, str) == True:
        low_temp = float(low_temp)
        print("Minimum temperature: " + str(low_temp))
    else:
        print("Error! Invalid minimum temperature!")
        exit()
		
    #print("Maximum temperture:")
    #high_temp = input()
	
    if isinstance(high_temp, str) == True:
        high_temp = float(high_temp)
        print("Maximum temperature: " + str(high_temp))
    else:
        print("Error! Invalid maximum temperature!")
        exit()
    
    header = []
    sequence = []
    no_header = []
    no_sequence = []

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #sequence = str(seq_record.seq).upper()
        if float(lines[i]) >= low_temp and float(lines[i]) <= high_temp:
            header.append(seq_record.description)
            sequence.append(seq_record.seq)
            #sequences[sequence] = seq_record.description
            i = i + 1
        else:
            no_header.append(seq_record.description)
            no_sequence.append(seq_record.seq)
            #no_seq[sequence] = seq_record.description
            #next
            i = i + 1

    # Create output file 

    # Create a file in the same directory where script was run
    with open("2_TmFiltered_" + fasta_file, "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for seq,head in zip(sequence,header):
            output_file.write(">" + str(head) + "\n" + str(seq) + "\n")
    with open("2_no_TmFiltered_" + fasta_file, "w+") as output_file:
        for seq,head in zip(no_sequence,no_header):
            output_file.write(">" + str(head) + "\n" + str(seq) + "\n")
	
    print("\n" + "Removed bad probes!\nOutput file: 2_TmFiltered_" + fasta_file + "\nUnwanted sequences found in output file: 2_no_TmFiltered_" + fasta_file)

melt_filter_range(userParameters, low_temp, high_temp)

