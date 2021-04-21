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



def delete_fastqs(sample, result_file_list, ini_section='delete_fastqs'):
	out_file=os.path.join("delete_fastqs", "done")
        #bam_lst = [os.path.join(os.path.dirname(result_file), "*bam*") for result_file in result_file_list] 
        #print(" ".join(bam_lst))
        #print(bam_lst)

#		eric_out=os.path.join("fusions", "ericscript", sample, "out"),
	return Job(
		result_file_list,
		[out_file],
		[],
		command="""\
rm -rf {fastq_folder} && rm -rf {tophat2_bam} && touch {out_file}""".format(
		fastq_folder=os.path.join("fusions", "gunzip_fastq"),
		tophat2_bam=os.path.join("fusions", "tophat2"),
		out_file=out_file
		),
		removable_files=[]
	)
#		fastq_folder=os.path.join("fusions", "gunzip_fastq", sample),
#		tophat2_bam=os.path.join("fusions", "tophat2", sample, "*.ba?"),

def delete_bams(result_file_list, topdir, ini_section='delete_fastqs'):                     
        out_file=os.path.join("delete_fastqs", "done")                                        
        #bam_lst = [os.path.join(os.path.dirname(result_file), "*bam*") for result_file in result_file_list] 
        #print(" ".join(bam_lst))
        #print(bam_lst)

#               eric_out=os.path.join("fusions", "ericscript", sample, "out"),                
        return Job(
                result_file_list,
                [out_file],
                [],
                command="""\
find {topdir}/fusions | grep -v metafusion | grep -v cff | grep 'bam\|sam' | xargs rm""".format(
                out_file=out_file,
                topdir=topdir
                ),
                removable_files=[]
        )
#find {topdir}/fusions | grep -v metafusion | grep -v cff | grep bam > {topdir}/bam_paths.txt""".format(
