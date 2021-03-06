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

def defuse(in1fastq, in2fastq, out_dir, keep_bam, config_file=None, ini_section='defuse'):

    other_options = config.param(ini_section, 'other_options', required=False)
    result_file = os.path.join(out_dir, "results.filtered.tsv")
    if keep_bam: keep_bam=str(1)
    else: keep_bam=str(0)
 
    return Job(
        [in1fastq, in2fastq, config_file if config_file else None],
        [result_file],
        [["defuse", "module_defuse"]],
        command="""\
OS_VERSION=$(cat /etc/centos-release) && if [[ $OS_VERSION == *"7"* ]]; then DEFUSE_CONFIG=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/defuse_config-CENTOS7.txt; else DEFUSE_CONFIG=/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/pipeline/config_reference_files/defuse_config.txt; fi &&\\
defuse.pl \\
  {other_options} \\
  -c $DEFUSE_CONFIG\\
  {in1fastq} \\
  {in2fastq} \\
  -o {out_dir} && \\
if [ {keep_bam} -eq 1 ];then ls -d {out_dir}/*|grep -v 'result\|bam' |xargs rm -rf; else ls -d {out_dir}/* |grep -v result |xargs rm -rf;fi; if ls core* 1> /dev/null 2>&1; then rm core*; else echo "no core dump files to remove"; fi""".format(
        other_options= other_options if other_options else "",
        in1fastq="  -1 " + in1fastq,
        in2fastq="  -2 " + in2fastq,
        keep_bam=keep_bam,
        out_dir= out_dir
        ),
        removable_files=[out_dir + "/reads.fqi", out_dir + "/reads.names", out_dir + "/reads.?.fastq"]

    )
#config_file= config_file if config_file else config.param(ini_section, 'defuse_config', required=True),
#ls -d {out_dir}/*|grep -v result|xargs rm -rf; if ls core* 1> /dev/null 2>&1; then rm core*; else echo "no core dump files to remove"; fi""".format(
