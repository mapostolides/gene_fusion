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


def filter_merged_summary_file(sum_file, num_captured_reads, ini_section='valfilter_cff_and_sample_enrichment', config_file=None):#, seq_len=None, ref_file=None, config_file=None, ini_section='repeat_filter'):
    other_options = config.param(ini_section, 'other_options', required=False)
    merged_summary_dir = os.path.join("fusion_reads_capture", "merged_summary_stats")
    filtered_sum_file = os.path.join(merged_summary_dir, "merged." + str(num_captured_reads) +  ".summary") 
    return Job(
        [sum_file],
        [filtered_sum_file],
        [['cff_conversion', 'module_fusiontools']],
        command="""\
awk -v cov="{num_captured_reads}" '$8 >= cov' {sum_file} | sort > {filtered_sum_file}""".format(
        sum_file=sum_file,
        filtered_sum_file=filtered_sum_file,
        num_captured_reads=str(num_captured_reads),
        ),
        removable_files=[],
        name="filter_merged_summary_file"
    )

#/hpf/largeprojects/ccmbio/mapostolides/reads_capture/stjude_reads_capture/filter_cff_file_using_validation_pipeline_output_SAMPLE_ENRICHMENT_A_TP_NT.py
def filter_cff_and_sample_enrichment(filtered_sum_file, cff_file, ini_section='repeat_filter'):
    num_captured_reads = config.param('valfilter_cff_and_sample_enrichment', 'num_captured_reads', type='int') 
    seq_len=config.param(ini_section, 'seq_len', type='int')

    ##cff_file_valfiltered = os.path.join(".".join([cff_file, "valfilter", str(num_captured_reads)]))
    # INPUT FILES
    cff_file_bwafilter = os.path.join(".".join([cff_file,"bwafilter", str(seq_len)]))
    # OUTPUT FILES
    cff_file_valfiltered = os.path.join(".".join([cff_file, "valfilter", str(num_captured_reads)]))
    cff_file_bwafilter_valfiltered = os.path.join(".".join([cff_file,"bwafilter", str(seq_len), "valfilter", str(num_captured_reads)]))

    return Job(
        [filtered_sum_file, cff_file, cff_file_bwafilter],
        [cff_file_valfiltered, cff_file_bwafilter_valfiltered],
        [['cff_conversion', 'module_fusiontools']],
        command="""\
filter_cff_and_sample_enrichment_using_reads_capture_output.py {filtered_sum_file} {cff_file} > {cff_file_valfiltered} &&
filter_cff_and_sample_enrichment_using_reads_capture_output.py {filtered_sum_file} {cff_file_bwafilter} > {cff_file_bwafilter_valfiltered}""".format(
        filtered_sum_file=filtered_sum_file,
        cff_file=cff_file,
        cff_file_bwafilter=cff_file_bwafilter,
        cff_file_valfiltered=cff_file_valfiltered,
        cff_file_bwafilter_valfiltered=cff_file_bwafilter_valfiltered
        ),
        removable_files=[],
        name="cff_file_valfilter_and_sample_enrichment"
    )

