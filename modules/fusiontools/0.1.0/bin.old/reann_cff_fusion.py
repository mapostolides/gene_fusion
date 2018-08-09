#!/usr/bin/env python
import sys
sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pygeneann
import sequtils
import pysam

cff_file = sys.argv[1]
#ref_fa = sys.argv[2]
ensbed = sys.argv[2]
# gene order annotation test
ref = sys.argv[3]
gene_ann = pygeneann.GeneAnnotation(ensbed)
	

n = 1	
for line in open(cff_file, "r"):
	fusion = pygeneann.CffFusion(line)
	fusion.ann_gene_order(gene_ann)
	#annotate fusion id and seq
	fusion.fusion_id = "F" + (str(n)).zfill(8)
	pygeneann.get_fusion_seq(fusion, ref, 100)

	fusion.check_codon(gene_ann, ref)
	print fusion.tostring()
	n += 1
