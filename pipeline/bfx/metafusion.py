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

def run_metafusion_singularity(out_dir_abspath, genome_fasta=None, gene_info=None, gene_bed=None, recurrent_bedpe=None, ini_section='metafusion'):
    num_tools = config.param('metafusion', 'num_tools', type='int')

    cff_dir_abspath = os.path.join(out_dir_abspath, "fusions", "cff")
    metafusion_outdir_abspath = os.path.join(out_dir_abspath, "fusions", "metafusion")

    cff = os.path.join("fusions", "cff", "merged.cff")
    code_dir="/hpf/largeprojects/ccmbio/mapostolides/MetaFusion"

    return Job(
        [cff],
        ["final.n" +str(num_tools) +".cluster"],
        [],
        command="""\
module load Singularity; \\
cd {code_dir}/scripts ;\\ 
singularity exec -B /home -B {code_dir} -B /tmp -B /localhd/tmp -B {cff_dir_abspath} -B {metafusion_outdir_abspath} {code_dir}/MetaFusion.simg \\
bash MetaFusion.sh --outdir {metafusion_outdir_abspath} \\
  --cff {cff} \\
  --gene_bed {gene_bed} \\
  --fusion_annotator \\
  --genome_fasta {genome_fasta} \\
  --gene_info {gene_info} \\
  --num_tools={num_tools} \\
  --recurrent_bedpe {recurrent_bedpe} \\
  --scripts {code_dir}/scripts""".format(
        cff=cff,
        code_dir=code_dir,
        metafusion_outdir_abspath=metafusion_outdir_abspath,
        cff_dir_abspath=cff_dir_abspath,
        num_tools=str(num_tools),
        gene_bed=gene_bed if gene_bed else config.param(ini_section, 'gene_bed', type='filepath'),
        genome_fasta=genome_fasta if genome_fasta else config.param(ini_section, 'genome_fasta', type='filepath'),
        gene_info=gene_info if gene_info else config.param(ini_section, 'gene_info', type='filepath'),
        recurrent_bedpe=recurrent_bedpe if recurrent_bedpe else config.param(ini_section, 'recurrent_bedpe', type='filepath')
        ),
        removable_files=[]
    )
