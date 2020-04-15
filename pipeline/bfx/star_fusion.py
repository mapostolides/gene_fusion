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

#def star_fusion(in1fastq, in2fastq, out_dir, CTAT_resource_lib, config_file=None, ini_section='star_fusion'):
def star_fusion(in1fastq, in2fastq, CTAT_resource_lib, out_dir, config_file=None, ini_section='star_fusion'):

    other_options = config.param(ini_section, 'other_options', required=False)
    result_file = os.path.join(out_dir, "star-fusion.fusion_predictions.abridged.tsv")
 
    return Job(
        [in1fastq, in2fastq, config_file if config_file else None],
        [result_file],
        [["star_fusion", "module_star_fusion"], ["star_fusion","module_perl"], ["star_fusion", "module_gcc"], ["star_fusion", "module_samtools"]],
        command="""\
STAR-Fusion --CPU 8 \\
  {other_options} \\
  --left_fq {in1fastq} \\
  --right_fq {in2fastq} \\
  --genome_lib_dir {CTAT_resource_lib} \\
  --output_dir {out_dir} """.format(
        other_options= other_options if other_options else "",
        in1fastq=in1fastq,
        in2fastq=in2fastq,
        CTAT_resource_lib=CTAT_resource_lib,
        out_dir= out_dir
        ),
        removable_files=[]

    )
#  --output_dir {out_dir} && ls -d {out_dir}/*|grep -v predictions | xargs rm -rf """.format(
#        command="""\
#STAR-Fusion --CPU 8 \\
#  {other_options} \\
#  --left_fq {in1fastq} \\
#  --right_fq {in2fastq} \\
#  --genome_lib_dir {CTAT_resource_lib} \\
#  --output_dir {out_dir} && \\
#ls -d {out_dir}/*|grep -v result|xargs rm -rf""".format(

#STAR-Fusion --CPU 8 \
#    --left_fq $left_fastq \
#    --right_fq $right_fastq \
#    --genome_lib_dir $CTAT_resource_lib \
#    --output_dir $output_dir

