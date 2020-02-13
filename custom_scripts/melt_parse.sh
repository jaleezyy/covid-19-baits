#!/bin/bash

# Run melt.pl and parse the data such that only Tm data outputted to a new file
# melt.pl run with given probe .fasta
# parsing occurs such that the data is first split between the hybrid software first run in melt.pl and the actual d* and Tm data
# The 4th column contains data for Tm which is separated and the header is removed
# An output file named with Tm_$SOURCE.txt is generated in the current directory

set -e # exit if pipeline returns non-zero status
set -u # unset variables == error
set -o pipefail # return value of last command to exit with non-zero status

# -s: SOURCE = .fasta file with probe sequences 

while getopts ":s:" option; do
	case "${option}" in
		s) SOURCE=$OPTARG;;
esac
done

melt.pl $SOURCE > "Tm_"$SOURCE".txt_intermediate"

count=$(grep -n "dG" "Tm_"$SOURCE".txt_intermediate" | awk -F ":" '{print $1}') # find line number where data splits between both halves of melt.pl
echo $count
total=$(wc -l "Tm_"$SOURCE".txt_intermediate" | awk -F " " '{print $1}') # get total number of lines so a range can be set
echo $total

# parse data to isolate Tm data
cat "Tm_"$SOURCE".txt_intermediate" | sed -n $count,$total'p' | awk '{print $4}' | tail -$(($count-1)) | head -$(($count-1)) > "Tm_"$SOURCE".txt"  

rm *"_intermediate"