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
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

#['run_star_seqr', 'module_starseqr_python'],
#            ['run_star_seqr', 'module_star']
def run(fastqs1, fastqs2, output_dir, sample_name, keep_bam):
    if not isinstance(fastqs1, list):
        fastqs1 = [fastqs1]
        
    if not isinstance(fastqs2, list):
        fastqs2 = [fastqs2]
    prefix="out"
    #output_file = os.path.join(output_dir + "_STAR-SEQR", sample_name + "_STAR-SEQR_candidates.txt")
    #output_dir = os.path.join("fusions", "star_seqr", sample.name)
    #output_file = os.path.join(output_dir, prefix + "_STAR-SEQR", prefix  + "_STAR-SEQR_candidates.txt")
    output_file = os.path.join(output_dir, prefix  + "_STAR-SEQR_candidates.txt")
    #fusions/star_seqr/smc_rna_sim45/out_STAR-SEQR/out_STAR-SEQR_candidates.txt

    if keep_bam: keep_bam=str(1)
    else: keep_bam=str(0)
    return Job(
        fastqs1 + fastqs2,
        [output_file],
        [
        ['run_star_seqr', 'module_star']
        ],
    #module purge;source /home/mapostolides/miniconda3/etc/profile.d/conda.sh; conda activate starseqr2;starseqr.py -t {threads} {options} \\
        command="""\
    module purge;source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh; conda activate starseqr;starseqr.py -t {threads} {options} \\
      -i {genome_build} \\
      -g {gene_annot} \\
      -r {reference} \\
      -1 {fastq1} \\
      -2 {fastq2} \\
      -p {output_dir}/{prefix} && ls -d {output_dir}/{prefix}_STAR-SEQR/* | grep -v 'candidates\|breakpoints\|bam' | xargs rm -rf;\\
      num=$(cat {output_dir}/out_STAR-SEQR.log  | grep "No candidates left to process" | wc -l ); if [ $num -gt 0 ];then  cp {output_dir}/{prefix}_STAR-SEQR/out_STAR-SEQR_breakpoints.txt {output_dir}/{prefix}_STAR-SEQR/out_STAR-SEQR_candidates.txt ;fi;  \\
       mv {output_dir}/{prefix}_STAR-SEQR/out_STAR-SEQR_candidates.txt {output_dir}/out_STAR-SEQR_candidates.txt""".format(
            genome_build=config.param('run_star_seqr', 'genome_build'),
            gene_annot=config.param('run_star_seqr', 'gene_annot'),
            reference=config.param('run_star_seqr', 'reference'),
            threads=config.param('run_star_seqr', 'threads', type='posint'),
            options=config.param('run_star_seqr', 'options'),
            fastq1=",".join(fastq1 for fastq1 in fastqs1),
            fastq2=",".join(fastq2 for fastq2 in fastqs2),
            output_dir=output_dir,
            prefix=prefix
        ),
    )
#      -p {output_dir}/{prefix} && ls -d {output_dir}/{prefix}_STAR-SEQR/* | grep -v 'candidates\|breakpoints\|bam' | xargs rm -rf;\\
#      -p {output_dir}/{prefix} && ls -d {output_dir}/{prefix}_STAR-SEQR/* | grep -v 'candidates\|breakpoints\|Chimeric.out.junction' | xargs rm -rf;\\
      #&& ls -d {output_dir}/{prefix}_STAR-SEQR | grep -v 'candidates\|breakpoints\|Chimeric.out.junction' | xargs rm -rf """.format(
#num=$(cat {output_dir}/out_STAR-SEQR.log  | grep "No candidates left to process" | wc -l ); if [ $num -gt 0 ];then  cp /hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/out_STAR-SEQR_candidates.header.txt {output_dir}/{prefix}_STAR-SEQR/out_STAR-SEQR_candidates.txt ;fi
