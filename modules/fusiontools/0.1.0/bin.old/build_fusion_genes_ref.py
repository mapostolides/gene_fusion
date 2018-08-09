#!/usr/bin/env python
import sys
sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pygeneann
import sequtils
import pysam

cff_file = sys.argv[1]
ensbed = sys.argv[2]
ref_fa = sys.argv[3]
# gene order annotation test
#ref = sys.argv[3]
gene_ann = pygeneann.GeneAnnotation(ensbed)
	
ref = pysam.FastaFile(ref_fa)
#fasta_f = open("tmp.fa", "w")
seq_dict ={}
for line in open(cff_file, "r"):
	fusion = pygeneann.CffFusion(line)

	gene1 = fusion.reann_gene1
	gene2 = fusion.reann_gene2
	fusion_id = fusion.fusion_id

	trans_seqs1 = pygeneann.build_transcript_sequences(gene_ann, fusion.reann_gene1, ref)	
	trans_seqs2 = pygeneann.build_transcript_sequences(gene_ann, fusion.reann_gene2, ref)	

	if not trans_seqs1 or not trans_seqs2:
		continue

	for trans_id in trans_seqs1:
		trans_seq = trans_seqs1[trans_id]
		print ">Gene1_" + gene1 + "_" + gene2 + "_" + trans_id + "_" + fusion_id
		print trans_seq
	for trans_id in trans_seqs2:
		trans_seq = trans_seqs2[trans_id]
		print ">Gene2_" + gene1 + "_" +  gene2 + "_" + trans_id + "_" + fusion_id
		print trans_seq
#fasta_f.close()
