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

def make_fusion_list(cff_file, tool, out_dir, config_file=None, ini_section='fusioninspector'):
        fusion_list_file = os.path.join(out_dir, tool + ".fusion_list.txt")
        fusiontools_dir = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts"
        return Job(
        [cff_file],
        [fusion_list_file],
        command="""\
  {fusiontools_dir}/generate_fusioninspector_fusions_list.sh \\
  {cff_file} {fusion_list_file}""".format(
        fusiontools_dir=fusiontools_dir,
        cff_file=cff_file,
        fusion_list_file=fusion_list_file 
        ),
        removable_files=[]

    )

#[["fusioninspector", "module_fusioninspector"]],
def fusioninspector(fastq1, fastq2, out_dir, fusion_list_file, config_file=None, ini_section='fusioninspector'):

    other_options = config.param(ini_section, 'other_options', required=False)
    result_file = os.path.join(out_dir, "finspector.fusion_predictions.final.abridged.FFPM")
    ctat_lib_dir="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play/ctat_genome_lib_build_dir"
    return Job(
        [fusion_list_file],
        [result_file],
        [
         ["fusioninspector", "module_samtools"], 
         ["fusioninspector","module_perl"], 
         ["fusioninspector", "module_bgzip"], 
         ["fusioninspector", "module_trinityseq"],
         ["fusioninspector", "module_gcc" ] #gcc required by STAR   
         ],
        command="""\
export PATH="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/STAR/STAR-2.6.1c/source:$PATH" && \\
/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/FusionInspector-run/FusionInspector/FusionInspector \\
  {other_options} \\
  --fusions {fusion_list_file} \\
  --genome_lib {ctat_lib_dir} \\
  --left_fq {fastq1} \\
  --right_fq {fastq1} \\
  --out_dir {out_dir} \\
  --out_prefix finspector \\
  --prep_for_IGV""".format(
        other_options= other_options if other_options else "",
        fusion_list_file=fusion_list_file,
        ctat_lib_dir=ctat_lib_dir,
        fastq1=fastq1,
        fastq2=fastq2,
        out_dir= out_dir
        ),
        #removable_files=[out_dir + "/reads.fqi", out_dir + "/reads.names", out_dir + "/reads.?.fastq"]
        removable_files=[]

    )


#filtered_out_dir = os.path.join("fusions", "cff_filtered", sample.name)
def filter_cff_calls_using_FI_results(FI_out_dir, filtered_out_dir, sample, tool, cff_dir):
    """
    Select fusions from .cff file that match FusionInspector output.
    """
    #EXISTING FILES/DIRECTORIES
    fusiontools_dir = "/hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts"
    renamed_cff_file = os.path.join(cff_dir, sample.name + "." + tool + ".cff.renamed") 
    FI_output_file = os.path.join(FI_out_dir, "finspector.fusion_predictions.final.abridged.FFPM")

    # TO BE CREATED
    FI_validated_cff_file = os.path.join(cff_dir, sample.name + "." + tool + ".cff.renamed.filtered")
    FI_unvalidated_cff_file = os.path.join(filtered_out_dir, sample.name + "." + tool + ".cff.renamed.unvalidated")
    FI_summary_file = os.path.join(filtered_out_dir, sample.name + "." + tool + ".FI_summary")

# python $scripts_path/filter_fusions_using_fusioninspector_output.py $FI_output_file $renamed_cff_file $FI_validated_cff_file $FI_unvalidated_cff_file $FI_summary_file $caller_name
    return Job(
        [input_cff_file, FI_output_file],
        [filtered_cff_file],
        command="""\
{fusiontools_dir}/filter_fusions_using_fusioninspector_output.py \\
  {FI_output_file} \\
  {renamed_cff_file} \\
  {FI_validated_cff_file} \\
  {FI_unvalidated_cff_file} \\
  {FI_summary_file}
  {caller_name} """.format(
        fusiontools_dir=fusiontools_dir,
        FI_output_file=FI_output_file,
        renamed_cff_file=renamed_cff_file,
        FI_validated_cff_file=FI_validated_cff_file,
        FI_unvalidated_cff_file=FI_unvalidated_cff_file,
        FI_summary_file=FI_summary_file,
        caller_name=tool
        ),
        removable_files=[]

    )


 
