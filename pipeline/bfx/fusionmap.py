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


def fusionmap(in1fastq, in2fastq, out_dir, top_dir, config_file=None, ini_section='fusionmap'):

    other_options = config.param(ini_section, 'other_options', required=False)
    out_prefix = "02_RNA"
    fastq_path = os.path.dirname(in1fastq)
    link1fastq = os.path.join(fastq_path, "tmplink_1.fastq")
    fastq_path = os.path.dirname(in2fastq)
    link2fastq = os.path.join(fastq_path, "tmplink_2.fastq")
    result_file = os.path.join(out_dir, "02_RNA.FusionReport.txt")
    #top_dir="/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modules"
# The below command can be used to create a new installation of FusionMap in the output directory
# For now I am trying to just use my local copy of FusionMap and not create a new one each run, as this seems to be creating conflicts
#cp -R -u -p /hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modules/FusionMap_2015-03-31 {top_dir}/ && \\
#cp -TRv /hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modules/FusionMap_2015-03-31 {top_dir}/  && \\
#/hpf/largeprojects/ccmbio/epigen/test_runs/run_fusion/output-July-2-2020/FusionMap_2015-03-31-TMP/

    return Job(
        [in1fastq, in2fastq, config_file if config_file else None],
        [result_file],
        [["fusionmap", "module_fusionmap"],
         ["DEFAULT", "module_mono"]],
        command="""\
cp -TRv /hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modules/FusionMap_2015-03-31 {out_dir}/FusionMap_2015-03-31  && \\
ln -sf $PWD/{in1fastq} {link1fastq} && \\
ln -sf $PWD/{in2fastq} {link2fastq} && \\
echo "\n<Files>\n{link1fastq}\n{link2fastq}\n<Output>\nTempPath={out_dir}/FusionMapTemp\nOutputPath={out_dir}\nOutputName={out_prefix}" \\
  >{out_dir}/tmp.cfg &&
cat {out_dir}/tmp.cfg \\
  {config_file} \\
  >{out_dir}/fusionmap.cfg && \\
/hpf/tools/centos6/mono/2.10.9/bin/mono \\
  /hpf/tools/centos6/Oshell/20170421/oshell.exe \\
  --semap {out_dir}/FusionMap_2015-03-31 \\
  Human.B37.3 \\
  RefGene \\
  {out_dir}/fusionmap.cfg \\
  {other_options} && \\
  rm -rf {out_dir}/FusionMap_2015-03-31""".format(
        other_options=" \\\n  " + other_options if other_options else "",
        config_file=config.param(ini_section, 'config_file'),
        in1fastq=in1fastq,
        in2fastq=in2fastq,
        link1fastq=link1fastq,
        link2fastq=link2fastq,
        top_dir=top_dir,
        out_dir=out_dir,
        out_prefix=out_prefix,
        fastq_path=fastq_path
        ),
        removable_files=[]
    )


# NOTES

# fusionmap application directory: /hpf/largeprojects/ccmbio/jiangyue/tools/FusionMap/FusionMap_2015-03-31, 
# was previously being used as first parameter to --semap flag
# instead, due to permission issues, a new instance of the FusionMap_2015-03-31 application is being generated in the top-level output directory of the pipeline
# each time it is run
