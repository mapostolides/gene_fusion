#!/usr/bin/env python
import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann
import sequtils
import pysam
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('cff_file', action='store', help='CFF file, can be .cff or cff.reann')
parser.add_argument('ensbed', action='store', help='Ensemble gene file')
parser.add_argument('ref_fa', action='store', help='Reference genome file')

args = parser.parse_args()

gene_ann = pygeneann.GeneAnnotation(args.ensbed)

n = 1	
for line in open(args.cff_file, "r"):
		
	fusion = pygeneann.CffFusion(line)
	orig_pos = fusion.pos1 # record the pos before shift
	fusion.ann_gene_order(gene_ann)
	#annotate fusion id and seq
	fusion.fusion_id = "F" + (str(n)).zfill(8)

	fusion.check_codon(gene_ann, args.ref_fa)
	if orig_pos != fusion.pos1: # the breakpoints have been shifted
		fusion.ann_gene_order(gene_ann) # remap genes with shifted breakpoints
	pygeneann.get_fusion_seq(fusion, args.ref_fa, 100)
	print fusion.tostring()
	n += 1
