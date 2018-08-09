#! /usr/bin/env python

import sys
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/pygeneann/pygenefusionann")
#sys.path.append("/hpf/largeprojects/ccmbio/jiangyue/Genap_ccm/pygenefusionann/")
import pygeneann
import argparse

# instantiate parser object
parser = argparse.ArgumentParser()

parser.add_argument('fusion_cluster_file', action='store', help='Fusion reannation file clustered by head/tail genes. (.reann.cluster file)')
parser.add_argument('sample_name_file', action='store', help='A file containing all sample names')

args = parser.parse_args()

class Sampleinfo():
    def __init__(self, line):
        self.sample_name, self.disease, self.sample_type = line.split()

def output_filtered_list(category_list):
    for c in category_list:
        #print c.gene1_on_bnd
        #print c.gene1_close_to_bnd
        #print c.gene2_on_bnd
        #print c.gene2_close_to_bnd
        c.out()
def output_cnt_table(fusion_list, group_tag="NA"):
    sample_type_list = ["Tumor", "Normal", "Tumor,Normal"]
    category_list = ["ReadThrough", "GeneFusion", "TruncatedCoding", "TruncatedNoncoding", "NoDriverGene", "SameGene"]
    print "Group:", group_tag
    print "Category\t" + "\t".join(sample_type_list)
    for category in category_list:
        print category ,
        for sample_type in sample_type_list:
            filtered_list = category_stats.filter_sample_type(fusion_list, sample_type)
            filtered_list = category_stats.filter_inferred_type(filtered_list, category)
            print len(filtered_list),
        print
def get_sample_type(sample):
    if sample.startswith("DIPG"):
        if sample.endswith("T"):
            sample_type = "Tumor"
        elif sample.endswith("N"):
            sample_type = "Normal"
    elif sample.startswith("SJH"):
            sample_type = "Tumor"
    elif sample.startswith("GTEX"):
        sample_type = "Normal"
    else:
        print >> sys.stderr, "Unkonwn sample:", sample
    return sample_type
def output_sample_fusion_cnt(fusion_list, sampleinfo_dict, group_tag):
    sample_list = []
    for fu in fusion_list:
        sample_list.extend(fu.samples)
    out_list = []
    for sample in sampleinfo_dict:
        out_list.append([sample, sampleinfo_dict[sample].sample_type, group_tag, sample_list.count(sample)])
    for li in sorted(out_list, key=lambda x:(x[1],x[0])):
        print "\t".join(map(str, li))

category_stats = pygeneann.CategoryFusionStats(args.fusion_cluster_file)

filtered_list = category_stats.category_list
#filtered_list = category_stats.filter_split_cnt(filtered_list, 3)

sampleinfo_dict = {}
for line in open(args.sample_name_file, 'r'):
    sampleinfo = Sampleinfo(line)
    sampleinfo_dict.setdefault(sampleinfo.sample_name, sampleinfo)
filter_sample_list = sampleinfo_dict.keys()
#filtered_list = category_stats.filter_tools_name(filtered_list, "stjude_method_Valid")
#filtered_list = category_stats.filter_samples(filtered_list, filter_sample_list)
#filtered_list = category_stats.filter_inferred_type(filtered_list, "ReadThrough")
#filtered_list = category_stats.filter_inferred_type(filtered_list, "GeneFusion")
filtered_list = category_stats.filter_tools_name(filtered_list, "ericscript")
print >> sys.stderr, "Total input category number:", len(filtered_list)
group = "Total"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)
#output_filtered_list(filtered_list)

#filtered_list = category_stats.filter_tools_num(filtered_list, 2)
#group = "2_tools"
#output_cnt_table(filtered_list, group)

#filtered_list = category_stats.filter_split_cnt(filtered_list, 5)
#filtered_list = category_stats.filter_span_cnt(filtered_list, 5)
#group = "split_span_5"
#output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, sampleinfo_dict, group)
#output_filtered_list(filtered_list)

"""
filtered_list = category_stats.filter_tools_num(filtered_list, 3)
group = "3_tools"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)

filtered_list = category_stats.filter_close_to_bnd(filtered_list)
group = "on_bnd"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)

filtered_list = category_stats.filter_dna_supp(filtered_list)
group = "dna_supp"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, filter_sample_list, group)
filtered_list = category_stats.filter_split_cnt(filtered_list, 20)
filtered_list = category_stats.filter_span_cnt(filtered_list, 20)
group = "split_span_20"
output_cnt_table(filtered_list, group)
#output_sample_fusion_cnt(filtered_list, args.sample_name_file, group)
"""

"""
Sample filters
"""
#filtered_list = category_stats.filter_recurrent(filter_list, 5)
#print filtered_list[0].disease
#filtered_list = category_stats.filter_disease(filtered_list, "DIPG")
#filtered_list = category_stats.filter_sample_type(filtered_list, "Tumor")
#filtered_list = category_stats.filter_sample_number(filtered_list, 3, "DIPG")
#filtered_list = category_stats.filter_split_cnt(filtered_list, 10)
#filtered_list = category_stats.filter_span_cnt(filtered_list, 10)
#filtered_list = category_stats.filter_tools_name(filtered_list, "defuse")
#filtered_list = category_stats.filter_tools_num(filtered_list, 2)
#filtered_list = category_stats.filter_dna_supp(filtered_list)
#filtered_list = category_stats.filter_inferred_type(filtered_list, "TruncatedCoding")
#filtered_list = category_stats.filter_close_to_bnd(filtered_list)
#filtered_list = category_stats.filter_recurrent(5)

#output_cnt_table(filtered_list)
#output_sample_fusion_cnt(filtered_list, args.sample_name_file)
#filtered_list = category_stats.filter_recurrent(filtered_list, 2)

#output_filtered_list(filtered_list)



