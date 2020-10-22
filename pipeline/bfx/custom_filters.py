#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

# generates statistics on the category/cluster file, generated my the merge_and_reannotate_cff_fusion step 

def covfilter(out_dir):

    seq_len = config.param('repeat_filter', 'seq_len', type='int')
    num_captured_reads = config.param('valfilter_cff_and_sample_enrichment', 'num_captured_reads', type='int')
    cov = config.param('fusion_stats', 'covfilter', type='int')

    # input files
    cluster_file = os.path.join(out_dir, "merged.cff.reann.dnasupp.cluster") 
    filter_bwa_cluster_file= os.path.join(out_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len) + ".cluster" )
    filter_val_cluster_file= os.path.join(out_dir, "merged.cff.reann.dnasupp" + ".valfilter." + str(num_captured_reads) + ".cluster" ) 
    filter_bwa_val_cluster_file= os.path.join(out_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len) + ".valfilter." + str(num_captured_reads) + ".cluster" )

    # output files
    cluster_file_cov = cluster_file + ".covfilter." + str(cov)
    cluster_file_bwa_cov = filter_bwa_cluster_file + ".covfilter." + str(cov)
    cluster_file_val_cov = filter_val_cluster_file + ".covfilter." + str(cov)
    cluster_file_bwa_val_cov = filter_bwa_val_cluster_file + ".covfilter." + str(cov)

#python $analysis_dir/category_fusion_stats_Covfilter5.py $cluster_file $sampleinfo > $cluster_file\.covfilter5
    return Job(
        [cluster_file, filter_bwa_cluster_file, filter_val_cluster_file, filter_bwa_val_cluster_file],
        [cluster_file_cov, cluster_file_bwa_cov, cluster_file_val_cov, cluster_file_bwa_val_cov],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
category_fusion_stats_Covfilter.py {cluster_file} {cov}  > {cluster_file_cov} && \\
category_fusion_stats_Covfilter.py {filter_bwa_cluster_file} {cov}  > {cluster_file_bwa_cov} && \\
category_fusion_stats_Covfilter.py {filter_val_cluster_file} {cov}  > {cluster_file_val_cov} && \\
category_fusion_stats_Covfilter.py {filter_bwa_val_cluster_file} {cov}  > {cluster_file_bwa_val_cov}""".format(
        cov=cov,
        cluster_file=cluster_file,
        filter_bwa_cluster_file=filter_bwa_cluster_file,
        filter_val_cluster_file=filter_val_cluster_file,
        filter_bwa_val_cluster_file=filter_bwa_val_cluster_file,
        cluster_file_cov=cluster_file_cov,
        cluster_file_bwa_cov=cluster_file_bwa_cov,
        cluster_file_val_cov=cluster_file_val_cov,
        cluster_file_bwa_val_cov=cluster_file_bwa_val_cov
        ),
        removable_files=[],
        name="covfilter." + str(cov)
    )


def custom_filters(cff_dir, cluster, ini_section='fusion_stats'):
    #config.param(ini_section, 'other_options', required=False)
    cff=os.path.join(cff_dir, "merged.cff.renamed.reann")
    blck_filtered_cluster = cluster + ".blck_filter"
    RT_caller2_filtered_cluster = blck_filtered_cluster + ".RT_filter.callerfilter2" 

    return Job(
        [cff, cluster],
        [RT_caller2_filtered_cluster, blck_filtered_cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator/blacklist_filter_recurrent_breakpoints.sh {cff} {cluster} {outdir}  \\
  > {blck_filtered_cluster} && \\
cat {blck_filtered_cluster} | grep -v ReadThrough | awk '$8 ~ /.,./' \\
  > {RT_caller2_filtered_cluster}\\
""".format(
        cff=cff,
        cluster=cluster,
        outdir=cff_dir,
        RT_caller2_filtered_cluster=RT_caller2_filtered_cluster,
        blck_filtered_cluster=blck_filtered_cluster
        ),
        removable_files=[]
    )
#cat {cluster} | grep -v ReadThrough | awk '$8 ~ /.,./' \\
#  > {RT_caller2_filtered_cluster} && \\
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator/blacklist_filter_recurrent_breakpoints.sh {cff} {RT_caller2_filtered_cluster} {outdir}  \\
#  > {blck_filtered_cluster}""".format(

def custom_filters_OLD(cff_dir, cluster, ini_section='fusion_stats'):
    #config.param(ini_section, 'other_options', required=False)
    cff=os.path.join(cff_dir, "merged.cff.renamed.reann")
    filtered_cluster = cluster + ".blck_filter.RT_filter.callerfilter2"

    return Job(
        [cff, cluster],
        [filtered_cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
/hpf/largeprojects/ccmbio/mapostolides/MODULES/FusionAnnotator/TEST_FusionAnnotator/blacklist_filter_recurrent_breakpoints.sh {cff} {cluster} {outdir} | grep -v ReadThrough | awk '$8 ~ /.,./' \\
  > {filtered_cluster}""".format(
        cff=cff,
        cluster=cluster,
        outdir=cff_dir,
        filtered_cluster=filtered_cluster
        ),
        removable_files=[]
    )




def filter_adjacent_noncoding(cff_dir, ini_section='fusion_stats'):
    #config.param(ini_section, 'other_options', required=False)
    filters=".blck_filter.RT_filter.callerfilter2"
    cluster = os.path.join(cff_dir, "merged.cff.renamed.reann.cluster" + filters) 
    filtered_cluster = cluster + ".ANC_filter"

    return Job(
        [cluster],
        [filtered_cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
filter_adjacent_noncoding.py {cluster} \\
  > {filtered_cluster}""".format(
        cluster=cluster,
        filtered_cluster=filtered_cluster
        ),
        removable_files=[]
    )

#python /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/output_ANC_RT_SG.py $cluster > $outdir/$(basename $cluster).ANC_RT_SG
def output_cis_sage_file(cff_dir, ini_section='fusion_stats'):
    #config.param(ini_section, 'other_options', required=False)
    cluster = os.path.join(cff_dir, "merged.cff.renamed.reann.cluster")
    cis_sage_cluster = cluster + ".ANC_RT_SG"

    return Job(
        [cluster],
        [cis_sage_cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
output_ANC_RT_SG.py {cluster} \\
  > {cis_sage_cluster}""".format(
        cluster=cluster,
        cis_sage_cluster=cis_sage_cluster 
        ),
        removable_files=[]
    )

def generate_category_count_table(cff_dir, out_dir, ini_section='fusion_stats' ):
    """
    Generates a tsv file with the gene fusion category counts, per category, per caller
    """
    cluster_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp.bwafilter.30.cluster")

    return Job(
        [cluster_file],
        [os.path.join(out_dir, "category_count_file.txt")],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
generate_category_table.py \\
  {cluster_file} \\
  {out_dir}""".format(
        cluster_file=cluster_file,
        out_dir=out_dir,
        ),
        removable_files=[]
    )

def generate_categories_barplot(fusion_stats_dir, ini_section='fusion_stats'):
    """
    Creates a barplot out of the "category_count_file.txt" file
    """
    category_count_file =  os.path.join(fusion_stats_dir, "category_count_file.txt")
    barplot_file = os.path.join(fusion_stats_dir, 'barplot_fusion_categories.html')
    
#generate_categories_barplot.py  $category_count_file $path 
    return Job(
        [category_count_file],
        [barplot_file],
        [["fusion_stats", "module_python_plotly"],["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
generate_categories_barplot.py \\
  {category_count_file} \\
  {output_dir}""".format(
        category_count_file=category_count_file,
        output_dir=fusion_stats_dir,
        ),
        removable_files=[]
    )




#SCRAP/NOTES

# category_fusion_stats.py testfiles/merged.cff.reann.dnasupp.bwafilter.30.cluster testfiles/sampleinfo

# cluster type, head gene, tail gene, max split read cnt, max spanning read cnt,  sample type, disease
#tools, inferred fusion category, head breakpoint on boundary, head breakpoint close to boundary, tail breakpoint on boundary,
#tail breakpoint close to boundary, dna support, samples

#1.  cluster type
#2.  head gene
#3.  tail gene
#4.  max split read cnt
#5.  max spanning read cnt
#6.  sample type
#7.  disease
#8.  tools
#9.  inferred fusion category
#10. head breakpoint on boundary
#11. head breakpoint close to boundary
#12. tail breakpoint on boundary
#13. tail breakpoint close to boundary
#14. dna support
#15. samples

