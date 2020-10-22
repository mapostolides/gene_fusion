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

def integrate(accepted_bam, unmapped_bam, out_dir, ini_section='integrate'):

    other_options = config.param(ini_section, 'other_options', required=False)
    breakpoints_file = os.path.join(out_dir, "breakpoints.tsv")
    reads_file = os.path.join(out_dir, "reads.txt")

    return Job(
        [accepted_bam, unmapped_bam],
        [breakpoints_file, reads_file],
        [["bwa", "module_bwa"]],
        command="""\
/hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/INTEGRATE-build/bin/Integrate fusion \\
  {other_options} \\
  /hpf/largeprojects/ccmbio/jiangyue/hg19_decoy/human_g1k_v37_decoy.fasta \\
  /hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/annot.ucsc.txt \\
  /hpf/largeprojects/ccmbio/jiangyue/hg19_decoy/integrate_index \\
  {accepted_bam} \\
  {unmapped_bam}""".format(
        other_options=" \\\n  " + other_options if other_options else "",
        accepted_bam=accepted_bam,
        unmapped_bam=unmapped_bam
        ),
        removable_files=[]
    )
#/hpf/largeprojects/ccmbio/mapostolides/MODULES/INTEGRATE/annot.ucsc.txt
#/hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/annot.ucsc.txt

#/hpf/largeprojects/ccmbio/mapostolides/MODULES/Integrate/INTEGRATE-build/bin/Integrate
#/hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/INTEGRATE-build/bin/Integrate fusion
def make_result_file(out_dir, ini_section='make_integrate_result_file'):

    other_options = config.param(ini_section, 'other_options', required=False)
    result_file = os.path.join(out_dir, "breakpoints.cov.tsv")
    reads_file = os.path.join(out_dir, "reads.txt")
    breakpoints_file=os.path.join(out_dir, "breakpoints.tsv")
    cov_file = os.path.join(out_dir, "cov.txt")

    return Job(
        [reads_file, breakpoints_file],
        [result_file],
        command="""\
{awk_cmd} && \\
paste {breakpoints_file} {cov_file} > {result_file}""".format(
        other_options=" \\\n  " + other_options if other_options else "",
        breakpoints_file=breakpoints_file,
        result_file=result_file,
        cov_file=cov_file,
        reads_file=reads_file,
        awk_cmd="""grep ^Fusion """ + reads_file + """|awk 'BEGIN{print "NUM_EN_RNA","NUM_SP_RNA"}{print $7,$9}' > """ + cov_file
        ),
        removable_files=[]
    )
"""
Integrate fusion (options) reference.fasta annotation.txt directory_to_bwt accepted_hits.bam unmapped.bam (dna.tumor.bam dna.normal.bam)

options: -cfn      integer : Cutoff of spanning RNA-Seq reads for fusions with non-canonical
                             exonic boundaries.                                                         default: 3
         -rt       float   : Normal dna / tumor dna ratio. If the ratio is less than
                             this value, then dna reads from the normal dna data set 
                             supporting a fusion candidates are ignored.                                default: 0.0
         -minIntra integer : If only having RNA reads, a chimera with two adjacent
                             genes in order is annotated as intra_chromosomal rather than 
                             read_through if the distance of the two genes is longer than
                             this value.                                                                default: 400000
         -minW     float   : Mininum weight for the encompassing rna reads on an edge.                  default: 2.0
         -mb       integer : See subcommand "mkbwt".
                             This value can be larger than used by mkbwt.                               default: 10000000
         -reads    string  : File to store all the reads.                                               default: reads.txt
         -sum      string  : File to store summary.                                                     default: summary.tsv
         -ex       string  : File to store exons for fusions with canonical exonic boundaries.          default: exons.tsv
         -bk       string  : File to store breakpoints                                                  default: breakpoints.tsv
         -bedpe    string  : File to store breakpoints in bedpe format                                  default: bk_fusion.bedpe
         -vcf      string  : File to store breakpoints in vcf format                                    default: bk_sv.vcf
         -bacc     integer : max difference between spanning reads and annotation to decide canonical.  default: 1
         -largeNum integer : if a gene shows greater or equal to this number, remove it from results.   default: 4
         -sample   string  : sample name                                                                default: sample



"""
