#!/usr/bin/env python2.7

import re

# CGTCTCAAGAT
# NNNNNNNNNNNNNNNNNNNN
# AATTTCTACTAAGTGTAGAT
# NNNNNNNNNNNNNNNNNNNN
# TTTTCGAGACG

left_adapter="TCCTTCAAGTGGTAAATGCGGCGTCTCAAGAT"
right_adapter="TTTTCGAGACGGTGTTGTCGTGACTGTCTATCCGT"
direct_repeat="AATTTCTACTAAGTGTAGAT"

counter = 1
for i in range(0,34):


	temp_file_name= r"targets_2August2017_3kpad_space"+str(i)+"_Single_gRNA_Stats.txt"

	dictionary = {}
	bsmbI_fw=r"CGTCTC"
	bsmbI_rev=r"GAGACG"

	with open(temp_file_name,'r') as f:

		for line in f:
			line_listed=line.strip("\n").split("\t")
			sample_name=line_listed[12]
			guide_seq=line_listed[5]
			MITscore=line_listed[18]
			# print guide_seq
			if not line_listed[0]=="gRNA_ID":
				if float(MITscore)>=90.0 and re.search(bsmbI_rev,guide_seq)==None and re.search(bsmbI_fw,guide_seq)==None:
					if sample_name not in dictionary:
						dictionary[sample_name] = {'left':[], 'right':[]}
						orientation = line_listed[11]
						dictionary[sample_name][orientation].append(line_listed)
					else:
						orientation = line_listed[11]
						dictionary[sample_name][orientation].append(line_listed)
				# else:
				# 	print [guide_seq,MITscore]

	for sample_name in dictionary:
		dictionary[sample_name]["left"].sort(key=lambda x: int(x[2]), reverse=True)
		dictionary[sample_name]["right"].sort(key=lambda x: int(x[2]), reverse=False)


	guidesleft_detailed = {}
	guidesright_detailed = {}

	for i in dictionary:

		guidesleft=[0]
		guidesright=[0]

		guidesleft_detailed[i] = []
		guidesright_detailed[i] = []

		counter_l = 0
		while len(guidesleft) < 6 and counter_l < len(dictionary[i]['left']):
			if abs(int(dictionary[i]['left'][counter_l][2]) - guidesleft[-1]) > 24:
				guidesleft.append(int(dictionary[i]['left'][counter_l][2]))
				guidesleft_detailed[i].append(dictionary[i]['left'][counter_l])
			counter_l += 1

		counter_r = 0
		while len(guidesright) < 6 and counter_r < len(dictionary[i]['right']):
			if abs(int(dictionary[i]['right'][counter_r][2]) - guidesright[-1]) > 24:
				guidesright.append(int(dictionary[i]['right'][counter_r][2]))
				guidesright_detailed[i].append(dictionary[i]['right'][counter_r])
			counter_r += 1

		if len(guidesright_detailed[i]) < 5 or len(guidesleft_detailed[i]) < 5:
			print i

	f1 = open('summary.txt', 'a')
	f2 = open('summary.bed', 'a')
	f3 = open('adaptors_added.txt', 'a')
	f4 = open('500_nc.txt', 'r')
	f5 = open("for_ordering.txt","a")
	for sample_name in guidesright_detailed:
		for left, right in zip(guidesleft_detailed[sample_name], guidesright_detailed[sample_name]):
			f1.write('\t'.join(map(str, left)) + '\t' + '\t'.join(map(str, right)) + '\t' + str(int(right[2]) - int(left[2])) + '\n')
			f2.write(str(left[1]) + '\t' + str(left[2]) + '\t' +str(left[3]) + '\n')
			f2.write(str(right[1]) + '\t' + str(right[2]) + '\t' +str(right[3]) + '\n')
			f3.write('>FASTA_LINE' + str(counter) + '\n')
			f5.write(left_adapter + str(left[5].strip('\n')) + direct_repeat + str(right[5].strip('\n')) + right_adapter + '\n')
			f3.write(left_adapter + str(left[5].strip('\n')) + direct_repeat + str(right[5].strip('\n')) + right_adapter + '\n')
			counter += 1
for lol in f4:
	g1 = lol.strip('\n').split()[0]
	g2 = lol.strip('\n').split()[1]
	f3.write('>FASTA_LINE' + str(counter) + '\n')
	f3.write(left_adapter + str(g1) + direct_repeat + str(g2) + right_adapter + '\n')
	f5.write(left_adapter + str(g1) + direct_repeat + str(g2) + right_adapter + '\n')
	counter += 1

