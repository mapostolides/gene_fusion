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

def merge_cff_fusion(input_cff_files, out_dir, annotation_file=None, reference_file=None, ini_section='merge_and_reannotate_cff_fusion'):

    merged_cff = os.path.join(out_dir, "merged.cff")

    return Job(
        input_cff_files,
        [merged_cff],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
cat {cff_files} > {out_dir}/merged.cff\\
  """.format(
        cff_files=" \\\n".join(input_cff_files),
        out_dir=out_dir
        ),
        removable_files=[]
    )

# generates cluster file, which is the final output file of the pipeline
def reannotate_cff_fusion(input_cff_files, out_dir, annotation_file=None, reference_file=None, ini_section='reannotate_cff_fusion'):

    merged_cff = os.path.join(out_dir, "merged.cff.renamed")

    return Job(
        input_cff_files,
        [merged_cff+".reann"],
        [["reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
reann_cff_fusion.py \\
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

def cluster_reann_dnasupp_file(outdir, ini_section='merge_and_reannotate_cff_fusion', repeat_filter_section='repeat_filter'):

    cff_prefix="merged.cff.renamed.reann"

    #num_captured_reads = config.param('valfilter_cff_and_sample_enrichment', 'num_captured_reads', type='int')
    #input files
    reann_cff = os.path.join(outdir, cff_prefix) 
    #output files
    cluster = reann_cff + ".cluster"

    return Job(
        [reann_cff],
        [cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
RUN_cluster_genes_breakpoints.sh {reann_cff} {outdir} > {cluster}""".format(
        reann_cff=reann_cff,
        outdir=outdir,
        cluster=cluster
        ),
        removable_files=[]
    )

