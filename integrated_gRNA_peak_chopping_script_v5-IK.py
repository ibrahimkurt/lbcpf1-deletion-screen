##### Cpf1 gRNA identifier/ranker software to make DNA excisions
# Input: Bed file, genome fasta file, various gRNA search parameters
# Output: Ranked gRNA pairs for each genomic site based on distance from "core"

##### Import libraries
print "Loading python dependencies and input genome and bed."

import re
import argparse
from bioutilities import Genome_2bit, Coordinate
import subprocess as sb
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

##### Time script
startTime = datetime.now()

##### Argument handeling
parser = argparse.ArgumentParser(description = 'Cpf1 gRNA identifier/ranker software to make DNA excisions. INPUTS: Genome, Bed File. | OUTPUTS: Ranked gRNA Pairs.')
parser.add_argument('-bed', help = 'Input bed file to make DNA excisions.')
parser.add_argument('-genome', help = 'Input genome 2bit file.')
parser.add_argument('-search_dist', default = 2500, help = 'Distance for gRNA search.')
parser.add_argument('-nuclease', default = 'cpf1', help = 'Specification of different CRISPR endonucleases (cas9 or cpf1).')
parser.add_argument('-weights', default = 'none', help = 'Custom weights for nuclease.')
parser.add_argument('-core', default = 'site_specific', help = 'Different core boundary calculations (site specific, mean, max).')
parser.add_argument('-out_dir', default = 'wd', help = 'Output directory specified by full path.')
parser.add_argument('-homopolymer', default = 3, help = 'Homopolymer limit in gRNA.')
parser.add_argument('-gc', default = [.25,.60], help = 'GC content range for gRNA.')
parser.add_argument('-pairs', default = 15, help = 'Number of unique gRNA pairs per site.')
parser.add_argument('-gRNA_score', default = 90, help = 'MIT gRNA score filter.')
parser.add_argument('-mm', default = 3, help = 'Number of mismatches to look for through CasOffinder.')
parser.add_argument('-mm0', default = 2, help = 'Number of off-target perfect matches to tolerate before suppressing MIT score to 0.')
parser.add_argument('-overlap', default = 'F', help = 'Whether the paired gRNAs can overlap in the unique paired gRNA calling.')
args = parser.parse_args()

##### Initialize arguments
bed_arg = args.bed
genome_arg = args.genome
search_dist = int(args.search_dist)
nuclease = args.nuclease
core_bound = args.core
homopolymer_lim = int(args.homopolymer)
gc_lb = float(args.gc[0])
gc_ub = float(args.gc[1])
pairs = int(args.pairs)
gRNA_score = float(args.gRNA_score)
custom_weights = args.weights
mm = int(args.mm)
mm0 = int(args.mm0)
overlap = args.overlap

##### Set output directory
if args.out_dir == 'wd':
	out_dir = str(os.getcwd())
else:
	if str(args.out_dir).endswith('/'):
		out_dir = str(args.out_dir[:-1])
	else:
		out_dir = str(args.out_dir)

##### Name output files
single_gRNA_out = out_dir + '/' + bed_arg.split('.')[0] + '_Single_gRNA_Stats.txt'
paired_gRNA_out = out_dir + '/' + bed_arg.split('.')[0] + '_Paired_gRNA_Stats.txt'
filtered_paired_gRNA_out = out_dir + '/' + bed_arg.split('.')[0] + '_Filtered_Paired_gRNA_Stats.txt'
casoffinder_in = out_dir + '/' + bed_arg.split('.')[0] + '_CasOFFinder_Input.txt'
casoffinder_out = out_dir + '/' + bed_arg.split('.')[0] + '_CasOFFinder_Output.txt'
failed_out = out_dir + '/' + bed_arg.split('.')[0] + '_Failed_Sites.txt'
total_gRNA_bed = out_dir + '/' + bed_arg.split('.')[0] + '_gRNAs_Total.bed'
filtered_gRNA_bed = out_dir + '/' + bed_arg.split('.')[0] + '_gRNAs_Filtered.bed'

##### Input genome
genome = Genome_2bit(genome_arg)

##### Input bed file
bed_input = open(bed_arg,'r')
bed = {}
peak_size_values = []
total_ids = []
for line in bed_input:
	line = line.strip('\n').split()
	chrom, start, stop = line[:3]
	if len(line[3:]) > 0:
		info = '_'.join(map(str, line[3:]))
		id_name = chrom + '_' + start + '_' + stop + '_' + info
	else:
		id_name = chrom + '_' + start + '_' + stop
	peak_size = int(stop) - int(start)
	peak_size_values.append(peak_size)
	if chrom not in bed:
		bed[chrom] = []
	bed[chrom].append([int(start), int(stop), id_name])
	total_ids.append(id_name)

##### Reverse complement function
def reverse_complement(sequence):
	sequence = sequence.upper()
	new_sequence = ''
	for base in sequence:
		if base == 'A':
			new_sequence += 'T'
		elif base == 'T':
			new_sequence += 'A'
		elif base == 'C':
			new_sequence += 'G'
		elif base == 'G':
			new_sequence += 'C'
	return(new_sequence[::-1])

##### Calculate GC content
def gc_content(sequence):
	sequence = sequence.upper()
	gc = float(sequence.count('G') + sequence.count('C'))
	total = float(len(sequence))
	gc_perc = gc/total
	return gc_perc

##### Identify homopolymers
def homopolymer(sequence):
	curr_val = 1
	curr_nuc = sequence[0]
	repeats = []
	for i in range(1,len(sequence)):
		if curr_nuc == sequence[i]:
			curr_val += 1
		else:
			repeats.append(curr_val)
			curr_val = 1
			curr_nuc = sequence[i]
	repeats.append(curr_val)
	return max(repeats)

##### PAM site initialization and position weights dependent on nuclease selection
if nuclease == 'cpf1':
	PAM_top = 'TTT[ACG]'
	PAM_bottom = '[TGC]AAA'
	gRNA_length = 24
	hitScoreM = [0.6221989685106383, 0.851063829787234, 0.8367170885106383, 0.8398649234042553, 0.8084330919148935, 0.8058788093617021, 0.7838360799999999, 0.5738456672340425, 0.5553549557446809, 0.851063829787234, 0.742949109787234, 0.7205955608510638, 0.851063829787234, 0.851063829787234, 0.851063829787234, 0.851063829787234, 0.7627932280851064, 0.7108900604255319, 0.35856887914893615, 0.0] #Cpf1 Scaled
	#hitScoreM = [0.731083788,1,0.983142579,0.986841285,0.949908883,0.946907601,0.921007394,0.674268659,0.652542073,1,0.872965204,0.846699784,1,1,1,1,0.896282043,0.835295821,0.421318433,0] #LbCpf1
elif nuclease == 'cas9':
	PAM_top = '[ATCG]GG'
	PAM_bottom = 'CC[ATCG]'
	gRNA_length = 23
	hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583] #Cas9

##### Option for custom weights
if custom_weights != 'none':
	hitScoreM = []
	with open(custom_weights, 'r') as w:
		for weight in w:
			weight = weight.strip('\n')
			hitScoreM.append(float(weight))
assert(len(hitScoreM)==20)

##### Initialize data storage for output
gRNA_l = {}
gRNA_r = {}
unreachable_sites = []
gRNA_cut_sites = [[],[],[]]
total_gRNA = open(total_gRNA_bed, 'w')

##### Establish "core" peak size
if core_bound == 'mean':
	core_size = int(float(sum(peak_size_values))/float(len(peak_size_values)))
elif core_bound == 'max':
	core_size = max(peak_size_values)
elif isinstance(core_bound, int):
	core_size = int(core_bound)

##### Perform gRNA search
print "Searching for gRNAs for each genomic site."

for chrom in bed:
	for site in bed[chrom]:
		start, stop, sample = site
		midpoint = int((int(stop) - int(start))/2.0) + int(start)
		peak_size = int((int(stop) - int(start))/2.0)
		if core_bound == 'site_specific':
			left_bound = int(start)
			right_bound = int(stop)
			core_size = int(right_bound) - int(left_bound)
		else:
			left_bound = midpoint - int(core_size/2.0)
			right_bound = midpoint + int(core_size/2.0)
		left = -peak_size + 1
		right = peak_size - 1
		for i in range(left, right + 1):
			gRNA_cut_sites[2].append(i)
		gRNA_l[sample], gRNA_r[sample] = [],[]
		flank_left = genome.extract_sequence(Coordinate(chrom, left_bound - search_dist, left_bound)).upper()
		flank_right = genome.extract_sequence(Coordinate(chrom, right_bound, right_bound + search_dist)).upper()
		seq_matches_ind_lt = [[m.start(), m.end()] for m in re.finditer(PAM_top, flank_left, re.IGNORECASE)]
		seq_matches_ind_lb = [[m.start(), m.end()] for m in re.finditer(PAM_bottom, flank_left, re.IGNORECASE)]
		seq_matches_ind_rt = [[m.start(), m.end()] for m in re.finditer(PAM_top, flank_right, re.IGNORECASE)]
		seq_matches_ind_rb = [[m.start(), m.end()] for m in re.finditer(PAM_bottom, flank_right, re.IGNORECASE)]
		if seq_matches_ind_lt:
			for match in seq_matches_ind_lt:
				start,stop = match
				if nuclease == 'cpf1':
					start_g = left_bound - (len(flank_left) - start)
					stop_g = left_bound - (len(flank_left) - (start + gRNA_length))
					gRNA_seq = flank_left[start:start + gRNA_length + 3].upper()
					cut_site = int(stop_g) - 1
					gRNA_cut_sites[0].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length + 3:
						entry = [chrom, start_g, stop_g, '+', gRNA_seq[4:gRNA_length], gRNA_seq[:4], cut_site, gRNA_seq[gRNA_length-1:], gc_content(gRNA_seq[4:24]), homopolymer(gRNA_seq[4:24]), 'left', sample]
						gRNA_l[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
				elif nuclease == 'cas9':
					start_g = left_bound - (len(flank_left) - (stop - gRNA_length))
					stop_g = left_bound - (len(flank_left) - stop)
					gRNA_seq = flank_left[stop - gRNA_length:stop].upper()
					cut_site = int(stop_g) - 6
					gRNA_cut_sites[0].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length:
						entry = [chrom, start_g, stop_g, '+', gRNA_seq[:gRNA_length-3], gRNA_seq[gRNA_length-3:gRNA_length],cut_site, 'NA', gc_content(gRNA_seq[:20]), homopolymer(gRNA_seq[:20]), 'left', sample]
						gRNA_l[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
		if seq_matches_ind_lb:
			for match in seq_matches_ind_lb:
				start,stop = match
				if nuclease == 'cpf1':
					start_g = left_bound - (len(flank_left) - (stop - gRNA_length))
					stop_g = left_bound - (len(flank_left) - stop)
					gRNA_seq = reverse_complement(flank_left[stop - gRNA_length - 3:stop])
					cut_site = int(start_g) + 1
					gRNA_cut_sites[0].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length + 3:
						entry = [chrom, start_g, stop_g, '-', gRNA_seq[4:gRNA_length], gRNA_seq[:4], cut_site, reverse_complement(gRNA_seq[gRNA_length-1:]), gc_content(gRNA_seq[4:24]), homopolymer(gRNA_seq[4:24]), 'left', sample]
						gRNA_l[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
				elif nuclease == 'cas9':
					start_g = left_bound - (len(flank_left) - start)
					stop_g = left_bound - (len(flank_left) - (start + gRNA_length))
					gRNA_seq = reverse_complement(flank_left[start:start + gRNA_length])
					cut_site = int(start_g) + 6
					gRNA_cut_sites[0].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length:
						entry = [chrom, start_g, stop_g, '-', gRNA_seq[:gRNA_length-3], gRNA_seq[gRNA_length-3:gRNA_length], cut_site, 'NA', gc_content(gRNA_seq[:20]), homopolymer(gRNA_seq[:20]), 'left', sample]
						gRNA_l[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
		if seq_matches_ind_rt:
			for match in seq_matches_ind_rt:
				start,stop = match
				if nuclease == 'cpf1':
					start_g = right_bound + start - 1
					stop_g = right_bound + (start + gRNA_length) - 1
					gRNA_seq = flank_right[start:start + gRNA_length + 3].upper()
					cut_site = int(stop_g) - 1
					gRNA_cut_sites[1].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length + 3:
						entry = [chrom, start_g, stop_g, '+', gRNA_seq[4:gRNA_length], gRNA_seq[:4], cut_site, gRNA_seq[gRNA_length-1:], gc_content(gRNA_seq[4:24]), homopolymer(gRNA_seq[4:24]), 'right', sample]
						gRNA_r[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
				elif nuclease == 'cas9':
					start_g = right_bound + (stop - gRNA_length) - 1
					stop_g = right_bound + stop - 1
					gRNA_seq = flank_right[stop - gRNA_length:stop].upper()
					cut_site = int(stop_g) - 6
					gRNA_cut_sites[1].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length:
						entry = [chrom, start_g, stop_g, '+', gRNA_seq[:gRNA_length-3], gRNA_seq[gRNA_length-3:gRNA_length], cut_site, 'NA', gc_content(gRNA_seq[:20]), homopolymer(gRNA_seq[:20]), 'right', sample]
						gRNA_r[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
		if seq_matches_ind_rb:
			for match in seq_matches_ind_rb:
				start,stop = match
				if nuclease == 'cpf1':
					start_g = right_bound + (stop - gRNA_length) - 1
					stop_g = right_bound + stop - 1
					gRNA_seq = reverse_complement(flank_right[stop - gRNA_length - 3:stop])
					cut_site = int(start_g) + 1
					gRNA_cut_sites[1].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length + 3:
						entry = [chrom, start_g, stop_g, '-', gRNA_seq[4:gRNA_length], gRNA_seq[:4], cut_site, reverse_complement(gRNA_seq[gRNA_length-1:]), gc_content(gRNA_seq[4:24]), homopolymer(gRNA_seq[4:24]), 'right', sample]
						gRNA_r[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')
				elif nuclease == 'cas9':
					start_g = right_bound + start - 1
					stop_g = right_bound + (start + gRNA_length) - 1
					gRNA_seq = reverse_complement(flank_right[start:start + gRNA_length])
					cut_site = int(start_g) + 6
					gRNA_cut_sites[1].append(cut_site - midpoint)
					if len(gRNA_seq) == gRNA_length:
						entry = [chrom, start_g, stop_g, '-', gRNA_seq[:gRNA_length-3], gRNA_seq[gRNA_length-3:gRNA_length], cut_site, 'NA', gc_content(gRNA_seq[:20]), homopolymer(gRNA_seq[:20]), 'right', sample]
						gRNA_r[sample].append(entry)
						total_gRNA.write(str(chrom) + '\t' + str(start_g) + '\t' + str(stop_g) + '\n')

total_gRNA.close()

##### Output single gRNAs stats
print "Outputting identified gRNAs for each genomic site."

total_gRNAs = []
single_out = open(single_gRNA_out, 'w')
counter = 1
single_out.write('\t'.join(['gRNA_ID','chr','start','stop','strand','gRNA_sequence','PAM', 'cut_site','sticky_end','GC','homopolymer','side','sample_type']) + '\n')
for chrom in gRNA_l:
	for gRNA in gRNA_l[chrom]:
		single_out.write(str(counter) + '\t' + '\t'.join(map(str, gRNA)) + '\n')
		counter += 1
		if gRNA[4] not in total_gRNAs:
			total_gRNAs.append(gRNA[4])
	for gRNA in gRNA_r[chrom]:
		single_out.write(str(counter) + '\t' + '\t'.join(map(str, gRNA)) + '\n')
		counter += 1
		if gRNA[4] not in total_gRNAs:
			total_gRNAs.append(gRNA[4])

single_out.close()

##### Create Cas-OFFinder input file
print "Performing off-target identification for all gRNAs."

co_in = open(casoffinder_in, 'w')
co_in.write(genome_arg + '\n')
if nuclease == 'cpf1':
	co_in.write('TTTVNNNNNNNNNNNNNNNNNNNN' + '\n')
	for gRNA in total_gRNAs:
		co_in.write('NNNN' + str(gRNA) + '\t' + str(mm) + '\n')
elif nuclease == 'cas9':
	co_in.write('NNNNNNNNNNNNNNNNNNNNNGG' + '\n')
	for gRNA in total_gRNAs:
		co_in.write(str(gRNA) + 'NNN' + '\t' + str(mm) + '\n')

co_in.close()

##### Run Cas-OFFinder on input file above
cmd = 'cas-offinder %s G %s' %(casoffinder_in, casoffinder_out)
sb.call(cmd, shell = True)

##### MIT Off-target scoring for LbCpf1

def calcMitScore(string1, string2, startPos=0):

	assert(len(string1)==len(string2)==20)
	dists = [] # distances between mismatches, for part 2
	mmCount = 0 # number of mismatches, for part 3
	lastMmPos = None # position of last mismatch, used to calculate distance
	score1 = 1.0
	for pos in range(0, len(string1)):
		if string1[pos]!=string2[pos]:
			mmCount+=1
			if lastMmPos!=None:
				dists.append(pos-lastMmPos)
			score1 *= 1-hitScoreM[pos]
			lastMmPos = pos
	# 2nd part of the score - distribution of mismatches
	if mmCount<2: # special case, not shown in the paper
		score2 = 1.0
	else:
		avgDist = sum(dists)/len(dists)
		score2 = 1.0 / (((19-avgDist)/19.0) * 4 + 1)
	# 3rd part of the score - mismatch penalty
	if mmCount==0: # special case, not shown in the paper
		score3 = 1.0
	else:
		score3 = 1.0 / (mmCount**2)
	score = score1 * score2 * score3 * 100
	return score

##### Parse off-targets by on-target
print "Parsing off-target output file."

OT_dict = {}
perf_match = {}
OT_count = {}
with open(casoffinder_out,'r') as OT:
	for line in OT:
		full = line.strip('\n').split()
		if nuclease == 'cpf1':
			gRNA_seq = full[0][4:]
			if gRNA_seq not in OT_dict:
				OT_dict[gRNA_seq] = [full[3][4:].upper()]
				OT_count[gRNA_seq] = [0]*(mm+2)
				OT_count[gRNA_seq][int(full[5])] += 1
			else:
				OT_dict[gRNA_seq].append(full[3][4:].upper())
				OT_count[gRNA_seq][int(full[5])] += 1
				OT_count[gRNA_seq][-1] += 1
		elif nuclease == 'cas9':
			gRNA_seq = full[0][:20]
			if gRNA_seq not in OT_dict:
				OT_dict[gRNA_seq] = [full[3][:20].upper()]
				OT_count[gRNA_seq] = [0]*(mm+2)
				OT_count[gRNA_seq][int(full[5])] += 1
			else:
				OT_dict[gRNA_seq].append(full[3][:20].upper())
				OT_count[gRNA_seq][int(full[5])] += 1
				OT_count[gRNA_seq][-1] += 1

##### Score individual gRNAs and output scores
print "Scoring all the gRNAs based on off-targets."

gRNA_scores = {}
total_scores = []
temp = open('temp', 'w')
single_gRNAs = open(single_gRNA_out, 'r')
next(single_gRNAs)

mm_labels = []
for i in range(0, mm+1):
	mm_labels.append('mm' + str(i))
mm_labels.append('total_mm')

temp.write('\t'.join(['gRNA_ID','chr','start','stop','strand','gRNA_sequence','PAM', 'cut_site','sticky_end','GC','homopolymer','side','sample_type']) + '\t' + '\t'.join(map(str, mm_labels)) + '\t' + 'MIT_agg_score' + '\n')

for line in single_gRNAs:
	line = line.strip('\n')
	full = line.split()
	gRNA_seq = full[5]
	gRNA_scores[full[0]] = []
	OT_count[gRNA_seq][0] = OT_count[gRNA_seq][0] - 1
	if int(OT_count[gRNA_seq][0]) <= mm0:
		if gRNA_seq in OT_dict:
			for ot in OT_dict[gRNA_seq]:
				single_ot_score = float(calcMitScore(gRNA_seq, ot))
				gRNA_scores[full[0]].append(single_ot_score)
		else:
			gRNA_scores[full[0]].append(0)
	else:
		gRNA_scores[full[0]].append(10**6)
	agg_score = 100/(100 + (sum(gRNA_scores[full[0]]) - 100))*100.0
	total_scores.append(agg_score)
	temp.write(line + '\t' + '\t'.join(map(str, OT_count[gRNA_seq])) + '\t' + str(agg_score) + '\n')

temp.close()
single_gRNAs.close()

sb.call('mv temp %s' %(single_gRNA_out), shell = True)

##### gRNA Pairing and filtering for unique pairs
print "Filtering paired gRNAs for each genomic site."

paired_gRNAs = open(paired_gRNA_out, 'w')
filtered_paired_gRNAs = open(filtered_paired_gRNA_out, 'w')
filtered_gRNA = open(filtered_gRNA_bed, 'w')
single_gRNAs = open(single_gRNA_out, 'r')
gRNA_side = {'left':{}, 'right':{}}
next(single_gRNAs)
for line in single_gRNAs:
	line = line.strip('\n')
	full = line.split()
	side = full[11]
	if full[12] not in gRNA_side[side]:
		gRNA_side[side][full[12]] = [line]
	else:
		gRNA_side[side][full[12]].append(line)

paired_gRNAs.write('\t'.join(['site','gRNA_L_ID','gRNA_L_seq','gRNA_L_gc','gRNA_L_homopolymer','gRNA_L_totalmm','gRNA_L_MITscore','gRNA_R_ID','gRNA_R_seq','gRNA_R_gc','gRNA_R_homopolymer','gRNA_R_totalmm','gRNA_R_MITscore','homology','del_length']) + '\n')
filtered_paired_gRNAs.write('\t'.join(['site','gRNA_L_ID','gRNA_L_seq','gRNA_L_gc','gRNA_L_homopolymer','gRNA_L_totalmm','gRNA_L_MITscore','gRNA_R_ID','gRNA_R_seq','gRNA_R_gc','gRNA_R_homopolymer','gRNA_R_totalmm','gRNA_R_MITscore','homology','del_length']) + '\n')

problem_sites = []
del_length_total = []
gRNA_pairs = {}

for site in total_ids:
	if site in gRNA_side['left']:
		for gRNAs_l in gRNA_side['left'][site]:
			gRNAs_l = gRNAs_l.strip('\n')
			full_l = gRNAs_l.split()
			id_l = full_l[0]
			chrom_l = str(full_l[1])
			start_l = str(full_l[2])
			stop_l = str(full_l[3])
			cut_site_l = int(full_l[7])
			sticky_end_l = full_l[8]
			gc_l = float(full_l[9])
			homopolymer_l = int(full_l[10])
			score_l = float(full_l[-1])
			if site in gRNA_side['right']:
				for gRNAs_r in gRNA_side['right'][site]:
					gRNAs_r = gRNAs_r.strip('\n')
					full_r = gRNAs_r.split()
					id_r = full_r[0]
					chrom_r = str(full_r[1])
					start_r = str(full_r[2])
					stop_r = str(full_r[3])
					cut_site_r = int(full_r[7])
					sticky_end_r = full_r[8]
					gc_r = float(full_r[9])
					homopolymer_r = int(full_r[10])
					score_r = float(full_r[-1])
					if sticky_end_l == sticky_end_r != 'NA':
						homology = 'True'
					else:
						homology = 'False'
					del_length = cut_site_r - cut_site_l
					paired_gRNAs.write('\t'.join(map(str, [site, id_l, full_l[5], gc_l, homopolymer_l, full_l[-2], score_l, id_r, full_r[5], gc_r, homopolymer_r, full_r[-2], score_r, homology, del_length])) + '\n')
					if (gc_lb <= gc_l <= gc_ub) and (gc_lb <= gc_r <= gc_ub) and (homopolymer_l <= homopolymer_lim) and (homopolymer_r <= homopolymer_lim) and (score_l >= gRNA_score) and (score_r >= gRNA_score):
						del_length_total.append(del_length)
						if site not in gRNA_pairs:
							gRNA_pairs[site] = [[site, id_l, full_l[5], gc_l, homopolymer_l, full_l[-2], score_l, id_r, full_r[5], gc_r, homopolymer_r, full_r[-2], score_r, homology, del_length, chrom_l, start_l, stop_l, chrom_r, start_r, stop_r]]
						else:
							gRNA_pairs[site].append([site, id_l, full_l[5], gc_l, homopolymer_l, full_l[-2], score_l, id_r, full_r[5], gc_r, homopolymer_r, full_r[-2], score_r, homology, del_length, chrom_l, start_l, stop_l, chrom_r, start_r, stop_r])

ids_used = []
overlap_check_l = [0]
overlap_check_r = [0]
for site in total_ids:
	pair_count = 1
	if site in gRNA_pairs:
		gRNA_pairs_sorted = sorted(gRNA_pairs[site], key = lambda x: (len(x[13]), x[14]))
		for gRNA_pair in gRNA_pairs_sorted:
			if int(pair_count) <= int(pairs):
				if overlap.upper() == 'T':
					counter += 1
					if (gRNA_pair[1] not in ids_used) and (gRNA_pair[7] not in ids_used):
						filtered_paired_gRNAs.write('\t'.join(map(str, gRNA_pair[:-6])) + '\n')
						filtered_gRNA.write('\t'.join(map(str, gRNA_pair[-6:-3])) + '\n')
						filtered_gRNA.write('\t'.join(map(str, gRNA_pair[-3:])) + '\n')
						ids_used.append(gRNA_pair[1])
						ids_used.append(gRNA_pair[7])
						pair_count += 1
				elif overlap.upper() == 'F':
					counter += 1
					if (gRNA_pair[1] not in ids_used) and (gRNA_pair[7] not in ids_used) and (abs(int(gRNA_pair[-5]) - overlap_check_l[-1]) > gRNA_length) and (abs(int(gRNA_pair[-2]) - overlap_check_r[-1]) > gRNA_length):
						filtered_paired_gRNAs.write('\t'.join(map(str, gRNA_pair[:-6])) + '\n')
						filtered_gRNA.write('\t'.join(map(str, gRNA_pair[-6:-3])) + '\n')
						filtered_gRNA.write('\t'.join(map(str, gRNA_pair[-3:])) + '\n')
						ids_used.append(gRNA_pair[1])
						ids_used.append(gRNA_pair[7])
						overlap_check_l.append(int(gRNA_pair[-5]))
						overlap_check_r.append(int(gRNA_pair[-2]))
						pair_count += 1
			else:
				break
		if int(pair_count) < 5:
			problem_sites.append([site, pair_count])
	else:
		problem_sites.append([site, 0])

paired_gRNAs.close()
filtered_paired_gRNAs.close()
filtered_gRNA.close()
single_gRNAs.close()

failed_sites = open(failed_out, 'w')
if len(problem_sites) > 0:
	for problem in problem_sites:
		failed_sites.write('\t'.join(map(str, problem)) + '\n')
failed_sites.close()

##### Output graphs
print "Plotting graphs for paired gRNA stats."

fig1 = sns.distplot(gRNA_cut_sites[0], kde_kws={"shade": True, "label":"Left CRISPR Cuts"})
fig1 = sns.distplot(gRNA_cut_sites[1], kde_kws={"shade": True, "label":"Right CRISPR Cuts"})
fig1 = sns.distplot(gRNA_cut_sites[2], kde_kws={"shade": True, "label":"Genomic Bed Regions"})
fig1.set(xlabel='Genome Index (Bed region centered at 0)', ylabel='Density', title = 'Distribution of gRNA Cut Sites Relative to Genome Targets')
fig1.figure.savefig(out_dir + '/' + bed_arg.split('.')[0] + 'gRNA_Cut_Sites.png', dpi = 1500)
fig1.figure.clf()

fig2 = sns.distplot(del_length_total, kde_kws={"shade": True})
fig2.set(xlabel='Deletion Length', ylabel='Density', title = 'Distribution of Deletion Lengths for Filtered gRNA Pairs')
fig2.figure.savefig(out_dir + '/' + bed_arg.split('.')[0] + 'Paired_gRNA_Deletion_Length_Distribution.png', dpi = 1500)
fig2.figure.clf()

fig3 = sns.distplot(total_scores)
fig3.set(xlabel='MIT gRNA Score', ylabel='Density', title = 'Distribution of MIT gRNA Scores', xlim=(0, 100))
fig3.figure.savefig(out_dir + '/' + bed_arg.split('.')[0] + 'gRNA_Score_Distribution.png', dpi = 1500)
fig3.figure.clf()

print datetime.now() - startTime
