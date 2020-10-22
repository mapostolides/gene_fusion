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

def run_metafusion_singularity(cff_dir, out_dir, genome_fasta=None, gene_info=None, gene_bed=None, recurrent_bedpe=None, ini_section='metafusion'):
    num_tools = config.param('metafusion', 'num_tools', type='int')
    cff = os.path.join(cff_dir, "merged.cff")
    metafusion_dir="/hpf/largeprojects/ccmbio/mapostolides/MetaFusion"

    return Job(
        [cff],
        ["final.n" +str(num_tools) +".cluster"],
        [],
        command="""\
module load Singularity; \\
cd {metafusion_dir}/scripts ;\\ 
singularity exec -B /home -B {metafusion_dir} -B /tmp -B /localhd/tmp -B {cff_dir} -B {out_dir} {metafusion_dir}/MetaFusion.simg \\
bash MetaFusion.sh --outdir {out_dir} \\
  --cff {cff} \\
  --gene_bed {gene_bed} \\
  --fusion_annotator \\
  --genome_fasta {genome_fasta} \\
  --gene_info {gene_info} \\
  --num_tools={num_tools} \\
  --recurrent_bedpe {recurrent_bedpe} \\
  --scripts {metafusion_dir}/scripts""".format(
        cff=cff,
        metafusion_dir=metafusion_dir,
        cff_dir=cff_dir,
        out_dir=out_dir,
        num_tools=str(num_tools),
        gene_bed=gene_bed if gene_bed else config.param(ini_section, 'gene_bed', type='filepath'),
        genome_fasta=genome_fasta if genome_fasta else config.param(ini_section, 'genome_fasta', type='filepath'),
        gene_info=gene_info if gene_info else config.param(ini_section, 'gene_info', type='filepath'),
        recurrent_bedpe=recurrent_bedpe if recurrent_bedpe else config.param(ini_section, 'recurrent_bedpe', type='filepath')
        ),
        removable_files=[]
    )


#topdir=/hpf/largeprojects/ccmbio/mapostolides
#fusiontools=$topdir/MetaFusion/scripts
##REFERENCE FILES FILES
##runs_dir=$topdir/MetaFusion/RUNS
#runs_dir=/hpf/largeprojects/ccmbio/mapostolides/Sandbox/RUNS
#mkdir $runs_dir
#gene_bed=$topdir/MetaFusion/reference_files/ens_known_genes.renamed.ENSG.bed
#gene_info=$topdir/MetaFusion/reference_files/Homo_sapiens.gene_info
#genome_fasta=$topdir/MetaFusion/reference_files/human_g1k_v37_decoy.fasta
#recurrent_bedpe=$topdir/MetaFusion/reference_files/blocklist_breakpoints.bedpe
#
#date=Oct-19-2020.MetaFusion_command
#echo trusight 
#outdir=$runs_dir/trusight.$date
#echo generating output in $outdir
#mkdir $outdir
#cff=$topdir/MetaFusion/test_data/trusight_cff/trusight.cff
#
#cd /hpf/largeprojects/ccmbio/mapostolides/MetaFusion/scripts
#singularity exec -B /home -B /hpf/largeprojects/ccmbio/mapostolides/MetaFusion -B /tmp -B /localhd/tmp -B /hpf/largeprojects/ccmbio/mapostolides/Sandbox/RUNS /hpf/largeprojects/ccmbio/mapostolides/MetaFusion/MetaFusion.simg bash MetaFusion.sh --outdir $outdir \
#                 --cff $cff  \
#                 --gene_bed $gene_bed \
#                 --fusion_annotator \
#                 --genome_fasta $genome_fasta \
#                 --gene_info $gene_info \
#                 --num_tools=2  \
#                 --recurrent_bedpe $recurrent_bedpe \
#                 --scripts $fusiontools
