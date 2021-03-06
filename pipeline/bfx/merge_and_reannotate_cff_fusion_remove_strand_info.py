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
# awk -F '\t' -v OFS='\t' { if ($11!="defuse") {$3="NA"; $6="NA"; print} }
def merge_cff_fusion(input_cff_files, out_dir, annotation_file=None, reference_file=None, ini_section='merge_and_reannotate_cff_fusion'):

    other_options = config.param(ini_section, 'other_options', required=False)
    merged_cff = os.path.join(out_dir, "merged.cff")
    # removes caller-provided strand sign information for all callers except defuse, replaces with NA
    #awk_contents="""{$3="NA"; $6="NA"; print}"""  
    awk_contents="""{ if ($11!="defuse") {$3="NA"; $6="NA"; print} }"""
    return Job(
        input_cff_files,
        [merged_cff],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
cat {cff_files} > {out_dir}/merged.cff-temp;\\
awk -F '\t' -v OFS='\t' '{awk_contents}' {out_dir}/merged.cff-temp > {out_dir}/merged.cff;\\
rm {out_dir}/merged.cff-temp\\
  """.format(
        cff_files=" \\\n".join(input_cff_files),
        out_dir=out_dir,
        awk_contents=awk_contents
        ),
        removable_files=[]
    )

# generates cluster file, which is the final output file of the pipeline
def reannotate_cff_fusion(input_cff_files, out_dir, annotation_file=None, reference_file=None, ini_section='reannotate_cff_fusion'):

    other_options = config.param(ini_section, 'other_options', required=False)
    merged_cff = os.path.join(out_dir, "merged.cff")

    return Job(
        input_cff_files,
        [merged_cff+".reann"],
        [["reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
reann_cff_fusion-TEST_desktop_sandbox.py \\
  {merged_cff} \\
  {annotation_file} \\
  {reference_file} \\
  > {merged_cff}.reann""".format(
        merged_cff=merged_cff,
        annotation_file=annotation_file if annotation_file else config.param(ini_section, 'annotation_file', type='filepath'),
        reference_file=reference_file if reference_file else config.param(ini_section, 'reference_file', type='filepath')
        ),
        removable_files=[]
    )

def cluster_reann_dnasupp_file(out_dir, ini_section='merge_and_reannotate_cff_fusion', repeat_filter_section='repeat_filter'):

    other_options = config.param(ini_section, 'other_options', required=False)
    reann_dnasupp_file = os.path.join(out_dir, "merged.cff.reann.dnasupp")

    # load seq_len used in repeat_filter step
    seq_len = config.param(repeat_filter_section, 'seq_len', type='int')
    repeat_filtered_file = os.path.join(out_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len))
    output_file = reann_dnasupp_file + ".cluster"
    repeat_filter_out_file=repeat_filtered_file + ".cluster"

    return Job(
        [reann_dnasupp_file, repeat_filtered_file],
        [output_file,repeat_filter_out_file],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
generate_common_fusion_stats.py {reann_dnasupp_file} > {out_file} && \\
generate_common_fusion_stats.py {repeat_filtered_file} > {repeat_filter_out_file}""".format(
        reann_dnasupp_file=reann_dnasupp_file,
        repeat_filtered_file=repeat_filtered_file,
        out_file=output_file,
        repeat_filter_out_file=repeat_filter_out_file
        ),
        removable_files=[]
    )

