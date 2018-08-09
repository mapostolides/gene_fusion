#!/usr/bin/env python

import sys
sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pysam
import sequtils
import os

def cigar_to_score(read_list):
	
	best_score = 0
	for read in read_list:
		total_score = 0
		if read.is_unmapped:
			continue
		'''		
		if read.cigartuples:
			for t in read.cigartuples:
				if t[0] == 0: # match
					score = 10
				elif 1 <= t[0] <= 2: # indel
					score = -1
				elif t[0] == 4: # softclipping
					score = -1
				else:
					score = 0
				total_score += score * int(t[1])
			if read.has_tag("NM"): # BWA alignment
				total_score -= read.get_tag("NM") * 10
			if read.has_tag("nM"): # STAR alignment
				total_score -= read.get_tag("nM") * 10
			#total_score -= read.get_tag("NM") * 10
		'''
		total_score = read.get_tag("AS")
		best_score = max(total_score, best_score)
	return best_score
				
def get_softclip_len(read):
	sc_len = 0
	for t in read.cigartuples:
		if t[0] == 4: # softclipping
			sc_len += int(t[1])
	return sc_len





captured_bam_file = sys.argv[1]
star_bam_file = sys.argv[2]
original_cram_file = sys.argv[3]
sample_info_file = sys.argv[4]

#110829_UNC11-SN627_0148_BD0DGAABXX.1.sorted.mdup.norRNA.cram    BRCA    TP      EXACT
sample_info_dict = {}
for line in open(sample_info_file, "r"):
	cram_name, disease, sample_type, sample, match_type = line.split()
	sample_info_dict.setdefault(cram_name, (sample, sample_type, disease))


#filename = os.path.basename(original_bam_file)
#tmp = filename.split(".")
cram = os.path.basename(original_cram_file) ## use bam file name as sample name
if cram in sample_info_dict:
	sample, sample_type, disease = sample_info_dict[cram]
else:
	print >> sys.stderr, sample, "not in sample_info file."
	sys.exit(1)

star_bam = pysam.AlignmentFile(star_bam_file, "rb")
captured_bam = pysam.AlignmentFile(captured_bam_file, "rb")

star_aln_dict = {}
for read in star_bam.fetch(until_eof=True):
	star_aln_dict.setdefault(read.query_name, []).append(read)
fusion_read_cnt_dict = {}
counted_read_dict = {} # captured reads can only be used once
for refname in captured_bam.references:
	#if refname != "ZNF43_MIR3916_F00000022_100":
	#if not "ITGB1" in refname:
	#	continue
	filterd_reads = []

	fusion_reads_dict = {} # reads from the same pair count once
	captured_reads_dict = {}
	tmp = refname.split("_")
	gene1 = tmp[0]
	gene2 = tmp[1]
	ref_type = tmp[2]
	bp = int(tmp[5])- 1
	outseqs = []


	print "------------------------------------------------------------------------------------------------------------"
	n_fusion = 0 # filtered fusion read cnt
	n_trans = 0 # filtered transcript read cnt
	n_fusion_captured = 0 # all captured fusion read cnt
	n_trans_captured = 0 # all captured fusion read cnt
	for read in captured_bam.fetch(refname):
		cigartuples = read.cigartuples
		clip_len_up = 0
		clip_len_down = 0
		if read.reference_start	< bp -10 and read.reference_end > bp + 10 and get_softclip_len(read) <= read.query_length * 0.05: # fusion support reads
			## score filter
			if read.is_secondary or read.is_supplementary:
				continue
			read_id = "_1" if read.is_read1 else "_2"
			star_hit = 0
			if read.query_name + read_id in star_aln_dict:
				star_hit = len(star_aln_dict[read.query_name + read_id])
				star_score = cigar_to_score(star_aln_dict[read.query_name + read_id])
			else:
				star_score = 0	
			captured_score = cigar_to_score([read])	
			if ref_type == "FUSION":
				n_fusion_captured += 1
			else:
				n_trans_captured += 1
			if captured_score < star_score:
				continue
			## output alignment
			if read.query_name in counted_read_dict:
				continue
			else:
				counted_read_dict.setdefault(read.query_name, "")
			fusion_reads_dict.setdefault(read.query_name, "")
			if ref_type == "FUSION":
				n_fusion += 1
			else:
				n_trans += 1
			if ref_type == "FUSION":	
				if cigartuples[0][0] == 4:
					clip_len_up = cigartuples[0][1]
				if cigartuples[-1][0] == 4:
					clip_len_down = cigartuples[-1][1]
				leading_space = " " * (30 - clip_len_up + read.reference_start)

				seq1 = read.query_sequence[:clip_len_up]
				seq2 = read.query_sequence[clip_len_up:bp-read.reference_start+clip_len_up]
				seq3 = read.query_sequence[bp-read.reference_start+clip_len_up:read.query_length-clip_len_down]
				seq4 = read.query_sequence[read.query_length-clip_len_down:]

				outseq = leading_space + " " +  seq1 + " " + seq2 + " " + seq3 + " " +  seq4
				outseqs.append((len(leading_space), outseq))
		
	#print "SUMMARY", sample, gene1, gene2, sample_type, "DIPG", len(fusion_reads_dict), len(captured_reads_dict)
	print "INFO", ref_type, sample, gene1, gene2, sample_type, disease, n_fusion, n_trans, n_fusion_captured, n_trans_captured

	key = gene1 + "_" + gene2
	fusion_read_cnt_dict.setdefault(key, []).append((n_fusion, n_trans, n_fusion_captured, n_trans_captured))

	for seq in sorted(outseqs, key=lambda x:x[0]):
		print seq[1]
for key in fusion_read_cnt_dict:
	gene1, gene2 = key.split("_")
	cnts = fusion_read_cnt_dict[key]
	sum_fusion =  sum([x[0] for x in cnts])
	sum_trans =  sum([x[1] for x in cnts])
	sum_fusion_captured = sum([x[2] for x in cnts])
	sum_trans_captured = sum([x[3] for x in cnts])

	print "SUMMARY", sample, disease, gene1, gene2, sample_type, sum_fusion, sum_trans, sum_fusion_captured, sum_trans_captured
	
star_bam.close()
captured_bam.close()







