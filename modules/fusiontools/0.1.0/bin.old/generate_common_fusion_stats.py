#!/usr/bin/env python
import sys
sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
import pygeneann
import sequtils

cff_file = sys.argv[1]
#ref_fa = sys.argv[2]
	

cffstats = pygeneann.CffFusionStats(cff_file)	
cffstats.generate_common_fusion_stats_by_genes(cff_file)
