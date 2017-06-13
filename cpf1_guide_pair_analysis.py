#!/usr/bin/env python2.7

from Bio import SeqIO
from Bio.Seq import Seq
import csv
from collections import OrderedDict
import sys
import argparse
import re

index_1 = "AATTTCTACTAAGTGTAGAT" #stores the index of the LbCpf1 direct repeat as the beginning of the casette
index_2 = "TTTTTTT" #stores the terminal T repeat as the end of the casette

patterned={}
unpatterned={}

def analyze_guide_pairs(input_file, fastq_file, output_file):
	"""
	creates a dictionary with guide counts from fastq_file, writes to output_file
	fastq_file: forward read fastq file
	output_file: csv file to write guide dictionary to
	dictionary: guide sequence as key, guide count as entry
	"""

	num_reads = 0 #total number of reads processed
	perfect_matches = 0 # guides with perfect match to library
	non_perfect_matches = 0 #number of guides without a perfect match to the library
	key_not_found = 0 #count of reads where key was not found

	# open library sequences and initiate dictionary of read counts for each guide
	try:
		with open(input_file, mode='rU') as infile:
			reader = csv.reader(infile)
			dictionary = {rows[0]:0 for rows in reader}
	except:
		print('could not open', input_file)

	# open fastq file
	try:
		handle = open(fastq_file, "rU")
	except:
		print("could not find fastq file")
		return

	readiter = SeqIO.parse(handle, "fastq")
	for record in readiter:
		num_reads += 1
		read_sequence = str.upper(str(record.seq))
		key_index_1 = read_sequence.find(index_1)
		key_index_2 = read_sequence.find(index_2)
		key_index_3 = read_sequence.rfind(index_1)
	

		if key_index_1 >= 0 and key_index_2 >=0 and ((key_index_2)-(key_index_1)>0) and ((key_index_3)-(key_index_1))==40:
			start_index = key_index_1 + len(index_1)
			checkpoint = read_sequence[(start_index + 60):(start_index + 67)]
			guide = read_sequence[start_index:(start_index + 60)]

			if checkpoint==index_2:
				if guide not in patterned:
					patterned[guide]=1
				elif guide in patterned:
					patterned[guide]+=1

		else:
			rev_comp_read_sequence = str(record.seq.reverse_complement())
			key_index_1 = rev_comp_read_sequence.find(index_1)
			key_index_2 = rev_comp_read_sequence.find(index_2)
			key_index_3 = rev_comp_read_sequence.rfind(index_1)

			if key_index_1 >= 0 and key_index_2 >=0 and ((key_index_2)-(key_index_1)>0) and ((key_index_3)-(key_index_1))==40:
				start_index = key_index_1 + len(index_1)
				checkpoint = rev_comp_read_sequence[(start_index + 60):(start_index + 67)]
				guide = rev_comp_read_sequence[start_index:(start_index + 60)]

				if checkpoint==index_2:
					if guide not in patterned:
						patterned[guide]=1
					elif guide in patterned:
						patterned[guide]+=1

			else:
				if read_sequence not in unpatterned:
					unpatterned[read_sequence]=1
				elif guide in unpatterned:
					unpatterned[read_sequence]+=1

	patterned_sorted = OrderedDict(sorted(patterned.items(), key=lambda t: t[0]))
	unpatterned_sorted = OrderedDict(sorted(unpatterned.items(), key=lambda t: t[0]))

	with open("patterned.seq", mode="w") as outfile1:
		a=1
		for guide in patterned_sorted:
			count = patterned_sorted[guide]
			outfile1.write(">"+str(a)+"\n")
			a+=1
			outfile1.write(guide+"\t"+"\n")

	with open("patterned.txt", mode="w") as outfile6:
			for guide in patterned_sorted:
				count = patterned_sorted[guide]
				outfile6.write(guide+"\t"+str(count)+"\n")

	with open("unpatterned.seq", mode="w") as outfile2:
		a=1
		for guide in unpatterned_sorted:
			count = unpatterned_sorted[guide]
			outfile2.write(">"+str(a)+"\n")
			a+=1
			outfile2.write(guide+"\t"+"\n")

	with open("unpatterned.txt", mode="w") as outfile7:
			for guide in unpatterned_sorted:
				count = unpatterned_sorted[guide]
				outfile7.write(guide+"\t"+str(count)+"\n")

	# create ordered dictionary with guides and respective counts and output as a csv file
	dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
	with open(output_file, 'w') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',')
		for guide in patterned_sorted:
			if guide in dict_sorted:
				count = patterned_sorted[guide]
				mywriter.writerow([guide,count])

	with open("unsorted_library_count.csv", "w") as csvfile2:
		mywriter=csv.writer(csvfile2, delimiter=",")
		with open(input_file, mode='rU') as infile2:
			reader = csv.reader(infile2)
			for row in reader:
				row2=row[0]
				if row2 in patterned_sorted:
					count = patterned_sorted[row2]
					mywriter.writerow([row2,count])
				elif row2 not in patterned_sorted:
					count = 0
					mywriter.writerow([row2,count])



	with open("minilibrary.seq", mode="w") as outfile3:
		a=1
		for guide in dict_sorted:
			outfile3.write(">"+str(a)+"\n")
			a+=1
			outfile3.write(guide+"\t"+"\n")

	with open("library_count.seq", mode="w") as outfile4:
		a=1
		for guide in patterned_sorted:
			if guide in dict_sorted:
				outfile4.write(">"+str(a)+"\n")
				a+=1
				outfile4.write(guide+"\t"+"\n")

	with open("unplanned_pairs.seq", mode="w") as outfile5:
		a=1
		for guide in patterned_sorted:
			if guide not in dict_sorted:
				outfile5.write(">"+str(a)+"\n")
				a+=1
				outfile5.write(guide+"\t"+"\n")

	with open("unplanned_pairs.txt", mode="w") as outfile8:
		for guide in patterned_sorted:
			if guide not in dict_sorted:
				count = patterned_sorted[guide]
				outfile8.write(guide+"\t"+str(count)+"\n")


	handle.close()
	return


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Analyze sequencing data for sgRNA library distribution')
	parser.add_argument('-f', '--fastq', type=str, dest='fastq_file',
						help='fastq file name', default='NGS.fastq')
	parser.add_argument('-o', '--output', type=str, dest='output_file',
						help='output file name', default='library_count.csv')
	parser.add_argument('-i', '--input', type=str, dest='input_file',
						help='input file name', default='library_sequences.csv')
	args = parser.parse_args()

	analyze_guide_pairs(args.input_file, args.fastq_file, args.output_file)
