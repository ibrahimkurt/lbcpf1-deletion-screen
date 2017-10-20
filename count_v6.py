#!/usr/bin/env python2.7

from Bio import SeqIO
from Bio.Seq import Seq
import csv
from collections import OrderedDict
import numpy as np
import sys
import argparse
import re

#KEY_REGION_START = 15 #start index of key region
#KEY_REGION_END = 45 #end index of key region
KEY_1 = "AATTTCTACTAAGTGTAGAT" #identifies sequence before guide to determine guide position
KEY_2 = "TTTTTTT"

#dict1={}
#dict2={}

patterned={}
unpatterned={}

def count_spacers(input_file, fastq_file, output_file, guide_g):
	"""
	creates a dictionary with guide counts from fastq_file, writes to output_file
	fastq_file: forward read fastq file
	output_file: csv file to write guide dictionary to
	dictionary: guide sequence as key, guide count as entry
	"""

	#prefix_for_files=fastq_file[]

	num_reads = 0 #total number of reads processed
	perfect_matches = 0 # guides with perfect match to library
	non_perfect_matches = 0 #number of guides without a perfect match to the library
	key_not_found = 0 #count of reads where key was not found

	# add 'G' to key sequence if included in library
	#if guide_g:
	#	global KEY
	#	KEY += "G"

	# open library sequences and initiate dictionary of read counts for each guide
	try:
		with open(input_file, mode='rU') as infile: #rU mode is necessary for excel!
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

	# process reads in fastq file
	readiter = SeqIO.parse(handle, "fastq")
	for record in readiter: #contains the seq and Qscore etc.
		num_reads += 1
		read_sequence = str.upper(str(record.seq))
		#print(read_sequence)
		#key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
		#print(key_region)
		key_index_1 = read_sequence.find(KEY_1)
		key_index_2 = read_sequence.find(KEY_2)
		key_index_3 = read_sequence.rfind(KEY_1)
		#print(key_index_fw)
		#print(key_index_rev)

		if key_index_1 >= 0 and key_index_2 >=0 and ((key_index_2)-(key_index_1)>0) and ((key_index_3)-(key_index_1))==40:
			#print((key_index_3)-(key_index_1))
			start_index = key_index_1 + len(KEY_1)
			checkpoint = read_sequence[(start_index + 60):(start_index + 67)]
			guide = read_sequence[start_index:(start_index + 60)]

			if checkpoint==KEY_2:
				if guide not in patterned:
					patterned[guide]=1
				elif guide in patterned:
					patterned[guide]+=1

		else:
			rev_comp_read_sequence = str(record.seq.reverse_complement())
			key_index_1 = rev_comp_read_sequence.find(KEY_1)
			key_index_2 = rev_comp_read_sequence.find(KEY_2)
			key_index_3 = rev_comp_read_sequence.rfind(KEY_1)

			if key_index_1 >= 0 and key_index_2 >=0 and ((key_index_2)-(key_index_1)>0) and ((key_index_3)-(key_index_1))==40:
				start_index = key_index_1 + len(KEY_1)
				checkpoint = rev_comp_read_sequence[(start_index + 60):(start_index + 67)]
				guide = rev_comp_read_sequence[start_index:(start_index + 60)]

				if checkpoint==KEY_2:
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

	with open(fastq_file+".patterned.seq", mode="w") as outfile1:
		a=1
		for guide in patterned_sorted:
			count = patterned_sorted[guide]
			outfile1.write(">"+str(a)+"\n")
			a+=1
			outfile1.write(guide+"\t"+"\n")

	with open(fastq_file+".patterned.txt", mode="w") as outfile6:
			for guide in patterned_sorted:
				count = patterned_sorted[guide]
				outfile6.write(guide+"\t"+str(count)+"\n")

	with open(fastq_file+".unpatterned.seq", mode="w") as outfile2:
		a=1
		for guide in unpatterned_sorted:
			count = unpatterned_sorted[guide]
			outfile2.write(">"+str(a)+"\n")
			a+=1
			outfile2.write(guide+"\t"+"\n")

	with open(fastq_file+".unpatterned.txt", mode="w") as outfile7:
			for guide in unpatterned_sorted:
				count = unpatterned_sorted[guide]
				outfile7.write(guide+"\t"+str(count)+"\n")

	# create ordered dictionary with guides and respective counts and output as a csv file
	dict_sorted = OrderedDict(sorted(dictionary.items(), key=lambda t: t[0]))
	with open(fastq_file+"."+output_file, 'wb') as csvfile:
		mywriter = csv.writer(csvfile, delimiter=',')
		for guide in patterned_sorted:
			if guide in dict_sorted:
				count = patterned_sorted[guide]
				mywriter.writerow([guide,count])

	with open(fastq_file+".unsorted_library_count.csv", "wb") as csvfile2:
		#x=0
		total=0
		mywriter=csv.writer(csvfile2, delimiter=",")
		with open(input_file, mode='rU') as infile2: #rU mode is necessary for excel!
			reader = csv.reader(infile2)
			for row in reader:
				#x+=1
				row2=row[0]
				#print(x)
				if row2 in patterned_sorted:
					count = patterned_sorted[row2]
					#print num_reads
					percent = 100*float(count)/float(num_reads)
					print percent
					mywriter.writerow([row2,count,percent])
				elif row2 not in patterned_sorted:
					count = 1
					percent = 100*float(count)/float(num_reads)
					mywriter.writerow([row2,count,percent])
				total+=percent
			print total	


	with open(fastq_file+".minilibrary.seq", mode="w") as outfile3:
		a=1
		for guide in dict_sorted:
			outfile3.write(">"+str(a)+"\n")
			a+=1
			outfile3.write(guide+"\t"+"\n")

	with open(fastq_file+".library_count.seq", mode="w") as outfile4:
		a=1
		for guide in patterned_sorted:
			if guide in dict_sorted:
				outfile4.write(">"+str(a)+"\n")
				a+=1
				outfile4.write(guide+"\t"+"\n")

	with open(fastq_file+".unplanned_pairs.seq", mode="w") as outfile5:
		a=1
		for guide in patterned_sorted:
			if guide not in dict_sorted:
				outfile5.write(">"+str(a)+"\n")
				a+=1
				outfile5.write(guide+"\t"+"\n")

	with open(fastq_file+".unplanned_pairs.txt", mode="w") as outfile8:
		for guide in patterned_sorted:
			if guide not in dict_sorted:
				count = patterned_sorted[guide]
				outfile8.write(guide+"\t"+str(count)+"\n")


	handle.close()

	#print num_reads

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
	parser.add_argument('-no-g', dest='guide_g', help='presence of guanine before spacer', action='store_false')
	parser.set_defaults(guide_g=True)
	args = parser.parse_args()

	count_spacers(args.input_file, args.fastq_file, args.output_file, args.guide_g)
