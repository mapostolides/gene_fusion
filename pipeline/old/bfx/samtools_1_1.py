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

# for convertion of our cram2.0 file to bam only
# "view" is the only command

#/hpf/tools/centos6/samtools/1.1/bin/samtools view -b $cram  -o $cram_bam
from core.config import *
from core.job import *


def view(input, output=None, ini_section='convert_cram_to_fastq'):
	other_options = config.param(ini_section, 'other_options', required=False)
	return Job(
		[input],
		[output],
		[],
		command="""\
export REF_PATH=/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Analysis/TCGA_CRAM_BAM/cram_ref && \\
/hpf/tools/centos6/samtools/1.1/bin/samtools view {other_options} \\
  {input}{output}""".format(
		other_options=other_options,
		input=input,
		output=" \\\n  > " + output if output else ""
		),
		removable_files=[output]
	)

