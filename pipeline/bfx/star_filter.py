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
def star_filter(star_bam, out_dir, ini_section='star_filter'):
	other_options = config.param(ini_section, 'other_options', required=False)
	out_fastq1 = os.path.join(out_dir, os.path.basename(star_bam) + ".star.filtered.pair1.fastq")
	out_fastq2 = os.path.join(out_dir, os.path.basename(star_bam) + ".star.filtered.pair2.fastq")
	return Job(
		[star_bam],
		[out_fastq1, out_fastq2],
		[["star_filter", "module_fusiontools"]],
		command="""\
filter_star_alignment.py {star_bam} {out_fastq1} {out_fastq2}""".format(
		star_bam=star_bam,	
		out_dir=out_dir,
		out_fastq1=out_fastq1,
		out_fastq2=out_fastq2
		)
	)
