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
def run(fastq1, fastq2, output_dir, ini_section='chimerascan'):
    index="/hpf/largeprojects/ccmbio/mapostolides/MODULES/run_chimerascan/chim_indexdir"
    return Job(
        [fastq1, fastq2],
        [], #OUTPUT FILE TO BE ADDED
        [
        #[]
        ],
        command="""\
    module purge;source /hpf/largeprojects/ccmbio/mapostolides/MODULES/miniconda3/etc/profile.d/conda.sh; conda activate chimerascan;chimerascan_run.py {index} {fastq1} {fastq2} {output_dir}""".format(
            index=index,
            fastq1=fastq1,
            fastq2=fastq2,
            output_dir=output_dir
        ),
    )
      #&& ls -d {output_dir}/{prefix}_STAR-SEQR | grep -v 'candidates\|breakpoints\|Chimeric.out.junction' | xargs rm -rf """.format(
