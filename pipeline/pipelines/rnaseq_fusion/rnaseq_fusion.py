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
import argparse
import collections
import logging
import os
import re
import sys
from collections import defaultdict

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

from bfx import samtools_1_1
from bfx import star_seqr 
from bfx import arriba
from bfx import star_fusion
from bfx import defuse
from bfx import fusionmap
from bfx import tophat2
from bfx import integrate
from bfx import ericscript
from bfx import gunzip
from bfx import merge_fastq
from bfx import cff_conversion
from bfx import rename_genes 
from bfx import fusioninspector_cluster as fusioninspector
from bfx import check_dna_support_before_next_exon
from bfx import merge_and_reannotate_cff_fusion
from bfx import repeat_filter
from bfx import custom_filters 
from bfx import fusion_stats
from bfx import validate_fusions
from bfx import delete_fastqs

import utils

from bfx import bedtools
from bfx import cufflinks
from bfx import differential_expression
from bfx import gq_seq_utils
from bfx import htseq
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import star
from bfx import bvatools
from bfx import rmarkdown
from pipelines import common
import utils

# READS CAPTURE IMPORTS
from bfx import filter_caputred_reads_DEV as filter_caputred_reads
from bfx import fusion_reads_capture_results_sum_DEV as fusion_reads_capture_results_sum
from bfx import build_fusion_and_head_gene_ref_DEV as build_fusion_and_head_gene_ref
from bfx import bwa_fusion_reads_capture_DEV as bwa_fusion_reads_capture
from bfx import extract_captured_reads_and_realn
from bfx import valfilter

log = logging.getLogger(__name__)

class RnaFusion(common.Illumina):
    """
    RNAFusion Pipeline
    ================
    The Gene Fusion pipeline identifies gene fusion events using RNA-seq FASTQ files. BAM/SAM files can also be used as input, which will then be converted to FASTQ format   
    
    Four separate tools detect fusion events: 
    [deFuse](https://sourceforge.net/p/defuse/wiki/DeFuse/), 
    [FusionMap](http://www.arrayserver.com/wiki/index.php?title=FusionMap), 
    [EricScript](https://sites.google.com/site/bioericscript/home), 
    and [INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home/).
   
    Note that the discord_read_trim parameter in the defuse configuration file should be set according to the mean fragment length.
    See [HERE](https://sourceforge.net/p/defuse/wiki/FAQ/) for more information. You can get the mean fragment length for your sample
    by running defuse, and then looking at the defuse log file in the job_output folder
 
    Tophat2 is used to generate precursor files for the INTEGRATE fusion detection tool.
     
    The fusion detection results are combined into one common file (.cff) that gives information about gene fusions including gene names, 
    type of fusion (ex. read through vs. gene fusion), and the tools that identified each fusion event. 
    Additionally, if DNA sequencing information is available for the samples of interest, 
    the Gene Fusion Pipeline can check for DNA support of gene fusions detected from RNA. 
    
    The RNAseq pipeline requires a sampleinfo file to be provided, which is a tab-delimited document with each sample 
    as a line with information about sample disease. The first column gives sample name, second column gives the disease name,
    third column tells whether the sample comes from a tumor (TP) or normal (NT) tissue. Additional columns can give further 
    information about samples.

    In addition, a dnabam file must be provided, which gives the name of .bam file(s) associated with RNA-seq sample.
    If there is no DNA sequencing associated with the sample, provide the name of an empty file    

    For validation, an optional file containing fusion gene pairs can be provided with the --valfile flag. Used
    only if the validate_fusions step is being run, and assumes that the input sequence data contains the fusions provided in 
    the validation file. This step tests the effectiveness of the pipeline in detecting known fusions.
    
    Notes:
    -integrate and fusionmap are the least computationally intensive, ericscript is more intensive, 
    and defuse is the most computationally intensive
    
    README IS INCOMPLETE BETWEEN THIS POINT...

    THIS SUMMARY SECTION IS NOT GENERATED BY THE FUSION PIPELINE. SHOULD A STEP BE ADDED?
    Finally, a summary html report is automatically generated by the pipeline at the end of the analysis.
    This report contains description
    of the sequencing experiment as well as a detailed presentation of the pipeline steps and results.
    Various Quality Control (QC) summary statistics are included in the report and additional QC analysis
    is accessible for download directly through the report. The report includes also the main references
    of the software tools and methods used during the analysis, together with the full list of parameters
    that have been passed to the pipeline main script.
    
    An example of the RNA-Seq report for an analysis on Public Corriel CEPH B-cell is available for illustration
    purpose only: [RNA-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/rnaReport/index.html).
    
    MORE INFORMATION ABOUT THE PIPELINE IS CURRENTLY NOT AVAILABLE. SHOULD THIS BE ADDED?
    [Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_RNA-Seq.pptx) is more
    [nfor]ation about the RNA-Seq pipeline that you may find interesting.

    ... AND THIS POINT.    
    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("--sampleinfo", help="sample info file", type=file)
        self.argparser.add_argument("--dnabam", help="DNA bam list", type=file)
        # add optional fusion validation file for pipeline validation mode
        self.argparser.add_argument("--valfile", required=False, help="fusion validation set file", type=file)
        #Class variables
        self.tool_list = ["star_seqr", "arriba", "star_fusion", "fusionmap", "ericscript", "integrate", "defuse"]
        super(RnaFusion, self).__init__()

    def star(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
        using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        project_index_directory = "reference.Merged"
        project_junction_file = os.path.join("alignment_1stPass", "AllSamples.SJ.out.tab")
        individual_junction_list=[]
        ######
        #pass 1 -alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_1stPass_directory = os.path.join("alignment_1stPass", readset.sample.name, readset.name)
            individual_junction_list.append(os.path.join(alignment_1stPass_directory,"SJ.out.tab"))

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_1stPass_directory,
                genome_index_folder=None,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else ""
            )
            job.name = "star_align.1." + readset.name
            jobs.append(job)
        
        ######
        jobs.append(concat_jobs([
        #pass 1 - contatenate junction
        star.concatenate_junction(
            input_junction_files_list=individual_junction_list,
            output_junction_file=project_junction_file
        ),
        #pass 1 - genome indexing
        star.index(
            genome_index_folder=project_index_directory,
            junction_file=project_junction_file
        )], name = "star_index.AllSamples"))

        ######
        #Pass 2 - alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_2ndPass_directory = os.path.join("alignment", readset.sample.name, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_2ndPass_directory,
                genome_index_folder=project_index_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else "",
                create_wiggle_track=True,
                search_chimeres=True,
                cuff_follow=True,
                sort_bam=True
            )
            job.input_files.append(os.path.join(project_index_directory, "SAindex"))
 
            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            # remove older symlink before otherwise it raise an error if the link already exist (in case of redo)
            if len(readset.sample.readsets) == 1:
                readset_bam = os.path.join(alignment_2ndPass_directory, "Aligned.sortedByCoord.out.bam")
                sample_bam = os.path.join("alignment", readset.sample.name ,readset.sample.name + ".sorted.bam")
                job = concat_jobs([
                    job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + os.path.relpath(readset_bam, os.path.dirname(sample_bam)) + " " + sample_bam, removable_files=[sample_bam])])

            job.name = "star_align.2." + readset.name
            jobs.append(job)

        report_file = os.path.join("report", "RnaSeq.star.md")
        jobs.append(
            Job(
                [os.path.join("alignment", readset.sample.name, readset.name, "Aligned.sortedByCoord.out.bam") for readset in self.readsets],
                [report_file],
                [['star', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable scientific_name="{scientific_name}" \\
  --variable assembly="{assembly}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    scientific_name=config.param('star', 'scientific_name'),
                    assembly=config.param('star', 'assembly'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="star_report")
        )

        return jobs

#    def run_star_fusion(self):
#        """
#        STAR-Fusion is a component of the Trinity Cancer Transcriptome Analysis Toolkit (CTAT). Based on the STAR aligner it identifies candidate fusion transcripts supported by Illumina reads. 
#        https://github.com/STAR-Fusion/STAR-Fusion/wiki
#        """
#
#        jobs = []
#        
#        left_fastqs = defaultdict(list)
#        right_fastqs = defaultdict(list)
#
#        for readset in self.readsets:
#            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + "-trimmed-")
#            
#            if  readset.run_type == "PAIRED_END":
#                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
#                if readset.fastq1 and readset.fastq2:
#                    candidate_input_files.append([readset.fastq1, readset.fastq2])
#                if readset.bam:
#                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
#                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
#            elif readset.run_type == "SINGLE_END":
#                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
#                if readset.fastq1:
#                    candidate_input_files.append([readset.fastq1])
#                if readset.bam:
#                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
#                [fastq1] = self.select_input_files(candidate_input_files)
#                fastq2 = None
#
#            else:
#                raise Exception("Error: run type \"" + readset.run_type +
#                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
#
#            left_fastqs[readset.sample.name].append(fastq1)
#            right_fastqs[readset.sample.name].append(fastq2)
#
#        for sample in self.samples:
#            output_dir = os.path.join("fusion", sample.name, "star_fusion")
#            
#            mkdir_job = Job(command="mkdir -p " + output_dir)
#  
#            job = concat_jobs([
#                mkdir_job,
#                star_fusion.run(left_fastqs[sample.name], right_fastqs[sample.name], output_dir)
#            ], name="run_star_fusion." + sample.name)
#            job.samples = [sample]
#            jobs.append(job)
#
#        return jobs

    def picard_sam_to_fastq(self):
        """
        Convert SAM/BAM files from the input readset file into FASTQ format
        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
        rerwritten from common.Illumina.picard_sam_to_fastq, make directory for this step under result folder in case the orginal bam file directory is not writtable
        """
        jobs = []
        for readset in self.readsets:
            # If readset FASTQ files are available, skip this step
            if not readset.fastq1:
                if readset.cram:
                    # convert cram to bam then to fastq. fastq and bam are saved on localhd
                    out_bam = os.path.join("$TMPDIR", os.path.basename(readset.cram)+".bam")
                    cram2bam_job = samtools_1_1.view(readset.cram, out_bam)
                    if readset.run_type == "PAIRED_END":
                        out_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                        fastq1 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", out_bam)))
                        fastq2 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", out_bam)))
                    else:
                        raise Exception("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                    picard_job = picard.sam_to_fastq(out_bam, fastq1, fastq2)
                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        cram2bam_job,
                        picard_job
                    ], name= "picard_sam_to_fastq." + readset.name)
                    jobs.append(job)
                elif readset.bam:
                    if readset.run_type == "PAIRED_END":
                        out_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                        fastq1 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)))
                        fastq2 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)))
                    else:
                        raise Exception("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                    picard_job = picard.sam_to_fastq(readset.bam, fastq1, fastq2)
                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        picard_job
                    ], name= "picard_sam_to_fastq." + readset.name)
                    jobs.append(job)
                else:
                    raise Exception("Error: BAM file not available for readset \"" + readset.name + "\"!")
        return jobs

    def gunzip_fastq(self):
        """
        Gunzip .fastq.gz files 
        """
        jobs = []
        for readset in self.readsets:
            out_dir = os.path.join("fusions", "gunzip_fastq", readset.sample.name)
            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = []
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    picard_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                    candidate_input_files.append([os.path.join(picard_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", readset.bam))), os.path.join(picard_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)))])
                if readset.cram:
                    picard_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                    candidate_input_files.append([os.path.join(picard_dir, os.path.basename(readset.cram)+".pair1.fastq.gz"), os.path.join(picard_dir, os.path.basename(readset.cram)+".pair2.fastq.gz")])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!")
            gunzip1_job = gunzip.gunzip_fastq(fastq1, out_dir)
            gunzip2_job = gunzip.gunzip_fastq(fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                gunzip1_job,
                gunzip2_job
            ], name="gunzip_fastq." + readset.sample.name + "." + readset.name)

            jobs.append(job)

        return jobs


    def merge_fastq(self):
        """
        Merge paired end fastqs of the same sample
        """
        jobs = []
        for sample in self.samples:
            if len(sample.readsets) > 1:
                input_dir = os.path.join("fusions", "gunzip_fastq", sample.name)
                fastq1_list = []
                fastq2_list = []
                for readset in sample.readsets:
                    if readset.bam:
                        fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq", readset.bam)))
                        fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq", readset.bam)))
                    if readset.fastq1:
                        if readset.fastq1.endswith(".gz"):
                            # input files are gzipped fastqs
                            fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq1)))
                            fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq2)))
                        else:
                            # input files are fastqs
                            fastq1 = os.path.join(input_dir, os.path.basename(readset.fastq1))
                            fastq2 = os.path.join(input_dir, os.path.basename(readset.fastq2))
                    fastq1_list.append(fastq1)
                    fastq2_list.append(fastq2)
                merge_fastq_job = merge_fastq.merge_fastq(fastq1_list, fastq2_list, input_dir)
                job = concat_jobs([
                    merge_fastq_job,
                ], name="merge_fastq." + sample.name)
                jobs.append(job)

        return jobs

    def select_input_fastq(self, sample):
        """
        Select input fastqs for fusion callers according to readset.
        This function is called in the gene fusion caller functions.
        """
        input_dir = os.path.join("fusions", "gunzip_fastq", sample.name)
          
        if len(sample.readsets) > 1: # sample has more than 1 readset, use merged fastq
            fastq1 = os.path.join(input_dir, "merged.pair1.fastq")
            fastq2 = os.path.join(input_dir, "merged.pair2.fastq")
        else:
            #input files are bams
            readset = sample.readsets[0]
            if readset.bam:
                fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq", readset.bam)))
                fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq", readset.bam)))
            if readset.fastq1:
                if readset.fastq1.endswith(".gz"):
                    # input files are gzipped fastqs
                    fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq1)))
                    fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq2)))
                else:
                    # input files are fastqs
                    fastq1 = os.path.join(input_dir, os.path.basename(readset.fastq1))
                    fastq2 = os.path.join(input_dir, os.path.basename(readset.fastq2))
            if readset.cram:
                fastq1 = os.path.join(input_dir, os.path.basename(readset.cram)+".pair1.fastq")
                fastq2 = os.path.join(input_dir, os.path.basename(readset.cram)+".pair2.fastq")
        print >> sys.stderr, fastq1
        print >> sys.stderr, fastq2
        return fastq1, fastq2

    def run_arriba(self):
        """
        """
    
        jobs = []

        left_fastqs = defaultdict(list)
        right_fastqs = defaultdict(list)
    
        for readset in self.readsets:
            trim_dir = os.path.abspath("trim")
            trim_file_prefix = os.path.join(trim_dir, readset.sample.name,readset.name + "-trimmed-")
    
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
    
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            left_fastqs[readset.sample.name].append(fastq1)
            right_fastqs[readset.sample.name].append(fastq2)
            
        for sample in self.samples:
            #fastq1, fastq2 = self.select_input_fastq(sample)
            #left_fastqs[sample.name].append(fastq1)
            #right_fastqs[sample.name].append(fastq2)
            output_dir = os.path.join("fusions", "arriba", sample.name)
        
            mkdir_job = Job(command="mkdir -p " + output_dir)
        
            chgdir_job = Job(command="cd " + output_dir)

            back_to_outdir_job = Job(command="cd " + self._output_dir)
        
            job = concat_jobs([
                mkdir_job,
                chgdir_job,
                arriba.run(left_fastqs[sample.name], right_fastqs[sample.name], self._output_dir, output_dir),
                back_to_outdir_job
            ], name="run_arriba." + sample.name)
            job.samples = [sample]
            jobs.append(job)
    
        return jobs

    def run_star_seqr(self):
        """
        RNA Fusion Detection and Quantification using STAR
        https://github.com/ExpressionAnalysis/STAR-SEQR
        """

        jobs = []
        for sample in self.samples:
            if len(sample.readsets) > 1: 
                raise Exception("Error: only one read set per sample allowed") 
            if sample.readsets[0].bam:#.bam input
                fastq_dir = os.path.join("fusions", "picard_sam_to_fastq", sample.name)     
                left_fastq = os.path.join(fastq_dir, sample.name + ".sorted.mdup.pair1.fastq.gz")
                right_fastq = os.path.join(fastq_dir, sample.name + ".sorted.mdup.pair2.fastq.gz")
            elif sample.readsets[0].fastq2 and sample.readsets[0].fastq2.split(".")[-1] == "gz":
                #print(sample.readsets[0].fastq2)
                #print(sample.readsets[0].fastq2.split(".")[-1])
                left_fastq = sample.readsets[0].fastq1 
                right_fastq = sample.readsets[0].fastq2 
            else:
                raise Exception("Error: only .bam and .fastq.gz inputs allowed")
            output_dir = os.path.join("fusions", "star_seqr", sample.name)

            job = concat_jobs([
                Job(command="mkdir -p " + output_dir),
                star_seqr.run(left_fastq, right_fastq, output_dir, sample.name)
            ], name="run_star_seqr." + sample.name)
        
            job.samples = [sample]
            jobs.append(job)

        return jobs


    def star_fusion(self):
        """
        Run STAR-Fusion to call gene fusions
        """
        jobs = []
        CTAT_resource_lib="/hpf/largeprojects/ccmbio/mapostolides/validate_fusion/test_star_star-fusion/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play/ctat_genome_lib_build_dir"
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "star_fusion", sample.name)
            #star_fusion_job = star_fusion.star_fusion(fastq1, fastq2, out_dir, CTAT_resource_lib)
            star_fusion_job = star_fusion.star_fusion(fastq1, fastq2, CTAT_resource_lib, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                star_fusion_job
            ], name="star_fusion." + sample.name)

            jobs.append(job)

        return jobs


    def defuse(self):
        """
        Run Defuse to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "defuse", sample.name)
            defuse_job = defuse.defuse(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                defuse_job
            ], name="defuse." + sample.name)

            jobs.append(job)

        return jobs


    def fusionmap(self):
        """
        Run FusionMap to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            # add pipeline top outpud dir as input to bfx fusionmap script
            # self._output_dir assigned from command line args in pipeline.py
            top_dir = self._output_dir

            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "fusionmap", sample.name)
            fusionmap_job = fusionmap.fusionmap(fastq1, fastq2, out_dir, top_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                fusionmap_job,
                Job(command="ls " + out_dir + "/02_RNA*")
            
            ], name="fusionmap." + sample.name)

            jobs.append(job)

        return jobs


    def ericscript(self):
        """
        Run EricScript to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "ericscript", sample.name)
            ericscript_job = ericscript.ericscript(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                Job(command="rm -r " + out_dir),
                ericscript_job
            ], name="ericscript." + sample.name)

            jobs.append(job)

        return jobs    
    
    def tophat2(self):
        """
        Run Tophat2 for Integrate. Determines accepted hits and unmapped reads, and outputs 
        corresponding .bam files required as input files for integrate step.
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join(self.output_dir, "fusions", "tophat2", sample.name)
            tophat2_job = tophat2.tophat2(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                tophat2_job
            ], name="tophat2." + sample.name)

            jobs.append(job)

        return jobs    

    def integrate(self):
        """
        Run Integrate to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            input_dir = os.path.join("fusions", "tophat2", sample.name)
            accepted_bam = os.path.join(self.output_dir, input_dir, "accepted_hits.bam")
            unmapped_bam = os.path.join(self.output_dir, input_dir, "unmapped.bam")

            out_dir = os.path.join("fusions", "integrate", sample.name)
            integrate_job = integrate.integrate(accepted_bam, unmapped_bam, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                Job(command="cd " + out_dir),
                integrate_job,
                Job(command="cd -")
            ], name="integrate." + sample.name)

            jobs.append(job)

        return jobs
    
    def integrate_make_result_file(self):
        """
        Merge infomation from breakpoints.tsv and reads.txt
        """
        jobs = []
        for sample in self.samples:
            input_dir = os.path.join("fusions", "integrate", sample.name)

            make_result_job = integrate.make_result_file(input_dir)
            job = concat_jobs([
                make_result_job
            ], name="integrate_make_result." + sample.name)

            jobs.append(job)

        return jobs    
        

    def convert_fusion_results_to_cff(self):
        """
        Convert fusion results of all 4 gene fusion callers to cff format
        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        job_list = [Job(command="mkdir -p " + out_dir)]
        sampleinfo_file = os.path.relpath(self.args.sampleinfo.name, self.output_dir)
        
        for sample in self.samples:
            
            #self.tool_list = ["star_seqr", "arriba", "star_fusion", "fusionmap", "ericscript", "integrate", "defuse"]
            star_seqr_result = os.path.join("fusions", "star_seqr", sample.name, "out_STAR-SEQR", "out_STAR-SEQR_candidates.txt")
            #print >> sys.stderr, star_seqr_result 
            arriba_result = os.path.join("fusions", "arriba", sample.name, "fusions.tsv")
            star_fusion_result = os.path.join("fusions", "star_fusion", sample.name, "star-fusion.fusion_predictions.abridged.tsv")
            defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
            fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
            ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
            integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")

            tool_results = [("star_seqr",star_seqr_result), ("arriba", arriba_result), ("star_fusion", star_fusion_result), ("defuse", defuse_result), ("fusionmap", fusionmap_result), ("ericscript", ericscript_result), ("integrate", integrate_result)]
            #tool_results = [("arriba", arriba_result), ("star_fusion", star_fusion_result), ("defuse", defuse_result), ("fusionmap", fusionmap_result), ("ericscript", ericscript_result), ("integrate", integrate_result)]
            #determine sample_type
            """
            sample_type = ""
            for contrast in self.contrasts:
                if sample in contrast.controls:
                    sample_type = "Normal"
                elif sample in contrast.treatments:
                    sample_type = "Tumor"
                if sample_type:
                    disease_name = contrast.name
                    break    
            if not sample_type:
                raise Exception("Error: sample " + sample.name + " not found in design file " + self.args.design.name)
            """   
            #convert caller output files to common fusion format(cff) 
            for tool, result_file in tool_results:
                job = cff_conversion.cff_convert(sample.name, result_file, sampleinfo_file, tool, out_dir)
                job.command = job.command.strip()
                job_list.append(job)
        job = concat_jobs(job_list, name="cff_conversion")
        jobs.append(job)
        return jobs


#    def fusioninspector(self):
#        """
#        Create fusion_list files for input to fusioninspector, and run fusioninspector on all samples for all callers separately 
#        """
#        jobs = []
#        cff_dir = os.path.join("fusions", "cff")
#        tool_list = ["star_seqr", "arriba", "star_fusion", "defuse", "fusionmap", "ericscript", "integrate"]
#        # CREATE FUSION LIST FILES 
#        cff_files = []
#        for sample in self.samples:
#            fusion_list_out_dir = os.path.join("fusions", "fusioninspector", "FI_fusion_list_files", sample.name)
#            job_list = []
#            job_list.append(Job(command="mkdir -p " + fusion_list_out_dir))
#            # create job_list for current sample
#            for tool in tool_list:
#                cff_file = os.path.join(cff_dir, sample.name + "." + tool + ".cff.renamed")
#                fusion_list_job = fusioninspector.make_fusion_list(cff_file, tool, fusion_list_out_dir)
#                job_list.append(fusion_list_job)
#            # concat jobs
#            job = concat_jobs(job_list, name= "make_fusion_list." + sample.name )
#            jobs.append(job)
#        # RUN FUSIONINSPECTOR
#        #select input fastq files 
#        for sample in self.samples:
#            fusion_list_out_dir = os.path.join("fusions", "fusioninspector", "FI_fusion_list_files", sample.name)
#            fastq1, fastq2 = self.select_input_fastq(sample)
#            # run fusioninspector for each tool separately
#            for tool in tool_list:
#                out_dir = os.path.join("fusions", "fusioninspector", "fusioninspector_output", sample.name, tool )
#                fusion_list_file = os.path.join(fusion_list_out_dir, tool + ".fusion_list.txt")
#                fusioninspector_job = fusioninspector.fusioninspector(fastq1, fastq2, out_dir, fusion_list_file)
#                job = concat_jobs([
#                    Job(command="mkdir -p " + out_dir),
#                    fusioninspector_job
#                ], name="fusioninspector." + sample.name + "." + tool)
#                jobs.append(job)
#        return jobs


    def fusioninspector(self):
        """
        Create fusion_list files for input to fusioninspector, and run fusioninspector on all samples for all callers separately 
        """
        jobs = []
        cff_dir = os.path.join("fusions", "cff")
        # CREATE FUSION LIST FILES 
        for sample in self.samples:
            fusion_list_out_dir = os.path.join("fusions", "fusioninspector", "FI_fusion_list_files", sample.name)
            job_list = []
            job_list.append(Job(command="mkdir -p " + fusion_list_out_dir))
            # create job_list for current sample
            cluster_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp.cluster")
            fusion_list_job = fusioninspector.make_fusion_list(cluster_file, sample.name, fusion_list_out_dir)
            job_list.append(fusion_list_job)
            # concat jobs
            job = concat_jobs(job_list, name= "make_fusion_list." + sample.name )
            jobs.append(job)
        # RUN FUSIONINSPECTOR
        #select input fastq files 
        for sample in self.samples:
            fusion_list_out_dir = os.path.join("fusions", "fusioninspector", "FI_fusion_list_files", sample.name)
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "fusioninspector", "fusioninspector_output", sample.name)
            fusion_list_file = os.path.join(fusion_list_out_dir, sample.name + ".fusion_list.txt")
            fusioninspector_job = fusioninspector.fusioninspector(fastq1, fastq2, out_dir, fusion_list_file)
            job = concat_jobs([
                    Job(command="mkdir -p " + out_dir),
                    fusioninspector_job
                ], name="fusioninspector." + sample.name)
            jobs.append(job)
        return jobs

    def filter_cff_calls_using_fusioninspector_results(self):
        jobs = []
        cff_dir = os.path.join("fusions", "cff")
        # FILTER CFF FUSIONS USING FUSIONINSPECTOR OUTPUT
        for sample in self.samples:
            job_list = []
            filtered_out_dir = os.path.join("fusions", "cff_filtered", sample.name)
            job_list.append([Job(command="mkdir -p " + filtered_outdir)])
            FI_out_dir = os.path.join("fusions", "fusioninspector", "fusioninspector_output", sample.name) 
            filter_job = fusioninspector.filter_cff_calls_using_FI_results(FI_out_dir, filtered_out_dir, sample, cff_dir)
            job_list.append(filter_job)
            job = concat_jobs(job_list, name= "filter_cff." + sample.name )
        return jobs

#    def filter_cff_calls_using_fusioninspector_results(self):
#        jobs = []
#        cff_dir = os.path.join("fusions", "cff")
#        # FILTER CFF FUSIONS USING FUSIONINSPECTOR OUTPUT
#        for sample in self.samples:
#            job_list = []
#            filtered_out_dir = os.path.join("fusions", "cff_filtered", sample.name)
#            job_list.append([Job(command="mkdir -p " + filtered_outdir)])
#            for tool in tool_list:
#                FI_out_dir = os.path.join("fusions", "fusioninspector", "fusioninspector_output", sample.name, tool ) 
#                filter_job = fusioninspector.filter_cff_calls_using_FI_results(FI_out_dir, filtered_out_dir, sample, tool, cff_dir)
#                job_list.append(filter_job)
#            job = concat_jobs(job_list, name= "filter_cff." + sample.name )
#        return jobs


    def merge_cff_fusion(self):
        """
        Merge all cff files into one single file
        """
        jobs = []
        cff_files = []
        cff_dir = os.path.join("fusions", "cff")
        out_dir = os.path.join("fusions", "cff")
        # put defuse .cff file last, which means inverted defuse calls will be always be "fusion2" in "generate_common_fusion_stats_by_breakpoints" function of pygeneann.py. This makes sense, since defuse is only one to make "flipped/inverted" calls. If defuse is not "fusion2" this results in errors in the case where defuse makes a flipped call
        #tool_list = ["star_seqr", "arriba", "star_fusion", "fusionmap", "ericscript", "integrate", "defuse"]
        for tool in self.tool_list:
            #cff_files.extend([os.path.join(cff_dir, sample.name + "." + tool + ".cff.renamed") for sample in self.samples])
            cff_files.extend([os.path.join(cff_dir, sample.name + "." + tool + ".cff") for sample in self.samples])
        merge_job = merge_and_reannotate_cff_fusion.merge_cff_fusion(cff_files, out_dir)
        
        job = concat_jobs([ merge_job ], name="merge_cff_fusion")
        jobs.append(job)
        return jobs

    def rename_genes(self):
        """ 
        Rename genes to consensus gene names using custom renaming script . This allows consistency in merging/categorizing downstream
        """
        jobs = []
        self.tool_list = ["star_seqr", "arriba", "star_fusion", "defuse", "fusionmap", "ericscript", "integrate"]
        out_dir = os.path.join("fusions", "cff")
        rename_genes_job = rename_genes.rename_cff_file_genes(os.path.join(out_dir, "merged.cff"), out_dir)
        rename_genes_job._name = "rename_cff_genes"
        #job = concat_jobs(rename_genes_job, name= "rename_cff_genes" )
        jobs.append(rename_genes_job) 

        return jobs

    
    def reannotate_cff_fusion(self):
        """
        Reannotate merged cff file with given annotation files
        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")

        merged_cff_file = os.path.join("fusions", "cff", "merged.cff") 
        reann_job = merge_and_reannotate_cff_fusion.reannotate_cff_fusion([merged_cff_file], out_dir)
        job = concat_jobs([reann_job], name="reannotate_cff_fusion")
        jobs.append(job)

        return jobs
    
        
    def check_dna_support_before_next_exon(self):
        """
        Check DNA support (pair clusters) until the start of next exon/utr
        """
        jobs = []
        dna_bam_list = os.path.abspath(self.args.dnabam.name)
        tmp_dir = os.path.join("fusions", "tmp")
        reann_file = os.path.join("fusions", "cff", "merged.cff.reann")
        dna_supp_job = check_dna_support_before_next_exon.check_dna_support_before_next_exon(reann_file, dna_bam_list, tmp_dir)
        job = concat_jobs([
            Job(command="mkdir -p " + tmp_dir),
            dna_supp_job    
        ], name="check_dna_support_before_next_exon")

        jobs.append(job)
        return jobs


    def repeat_filter(self):
        """
        Filter fusions with repetitive boundary sequences by realigning a certain length of sequnces with BWA
        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        cff = os.path.join(out_dir, "merged.cff.reann.dnasupp")
        job = repeat_filter.repeat_filter(cff, out_dir)
        
        job = concat_jobs([
            job    
        ], name="repeat_filter")

        jobs.append(job)
        return jobs

##START fusion reads capture pipeline
    def build_fusion_and_head_gene_ref(self):
        """
        Build fusion reference together with the head gene's all transcripts sequences
        """

        jobs = []
        cff_dir = os.path.join("fusions", "cff")
        #cff_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp.bwafilter.30")
        cff_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp")
        out_dir = os.path.join("fusion_reads_capture", "fusion_refs")

        build_job = build_fusion_and_head_gene_ref.build_fusion_and_head_gene_ref(cff_file, out_dir)
        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            build_job
        ], name="build_fusion_and_head_gene_ref")

        jobs.append(job)
        return jobs

    def fastq_conversion_and_reads_capture(self):
        """
        Coert cram2.0 file to fastq file with samtools1.1 and picard
        """
        jobs = []
        # make .cff file of repeat_filter step input to this step
        seq_len = config.param('repeat_filter', 'seq_len', type='int')
        cff_dir = os.path.join("fusions", "cff")
        #cff_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len))
        cff_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp")
        ref = os.path.join("fusion_reads_capture", "fusion_refs", os.path.basename(cff_file)+".fa")
        #TODO make sure cram file works
        #cram_file = self.args.cff.name
        #out_dir = os.path.join("fusion_reads_capture", "cram_fastq")
        for readset in self.readsets:
            out_dir = os.path.join("fusion_reads_capture", "captured_bam", readset.sample.name)
            if not readset.fastq1:
                if readset.bam:
                    # convert bam to fastq, fastq saved on localhd
                    fastq1 = os.path.join("$TMPDIR", os.path.basename(readset.bam) + ".1.fastq")
                    fastq2 = os.path.join("$TMPDIR", os.path.basename(readset.bam) + ".2.fastq")
                    bam2fastq_job = picard.sam_to_fastq(readset.bam, fastq1, fastq2)
                    # bwa aln fastqs to capture reference
                    out_bam = os.path.join(out_dir, "captured.bam")
                    capture_job = bwa_fusion_reads_capture.bwa_fusion_reads_capture(fastq1, fastq2, ref, out_bam, read_group=None, ini_section='bwa_fusion_reads_capture')

                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        bam2fastq_job,
                        capture_job
                    ], name="fastq_conversion_and_reads_capture."+readset.sample.name)
                    #], name="fastq_conversion_and_reads_capture" )
                                        # manually set I/O for job
                    job._input_files = [readset.bam, ref]
                    job._output_files = [out_bam]
                    jobs.append(job)

                elif readset.cram:
                    # convert cram to bam then to fastq, fastq and bam are saved on localhd
                    out_bam = os.path.join("$TMPDIR", os.path.basename(readset.cram)+".bam")
                    fastq1 = out_bam + ".1.fastq"
                    fastq2 = out_bam + ".2.fastq"
                    #cram2bam_job = samtools_1_1.view(readset.cram, out_bam, "convert_cram_to_fastq")
                    cram2bam_job = samtools_1_1.view(readset.cram, out_bam)
                    bam2fastq_job = picard.sam_to_fastq(out_bam, fastq1, fastq2)

                    # bwa aln fastqs to capture reference
                    out_bam = os.path.join(out_dir, "captured.bam")                                                          
                    capture_job = bwa_fusion_reads_capture.bwa_fusion_reads_capture(fastq1, fastq2, ref, out_bam, read_group=None, ini_section='bwa_fusion_reads_capture')

                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        cram2bam_job,
                        bam2fastq_job,
                        capture_job
                    ], name="fastq_conversion_and_reads_capture" +readset.sample.name)
                    #], name="fastq_conversion_and_reads_capture"+readset.sample.name)
                                        # manually set I/O for job
                    job._input_files = [readset.cram, ref]
                    job._output_files = [out_bam]

                    jobs.append(job)
                else:
                    raise Exception("Error: CRAM file not available for readset \"" + readset.name + "\"!")
            else:
                fastq1 = readset.fastq1
                fastq2 = readset.fastq2

                # bwa aln fastqs to capture reference
                out_bam = os.path.join(out_dir, "captured.bam")
                capture_job = bwa_fusion_reads_capture.bwa_fusion_reads_capture(fastq1, fastq2, ref, out_bam, read_group=None, ini_section='bwa_fusion_reads_capture')

                job = concat_jobs([
                    Job(command="mkdir -p " + out_dir),
                    capture_job
                ], name="bwa_fusion_reads_capture."+readset.sample.name)

                jobs.append(job)

        return jobs

    def extract_captured_reads_and_realn(self):
        """
        BWA mem realign captured reads to hg + transcript junction references
        """

        jobs = []
        for readset in self.readsets:

            captured_bam = os.path.join("fusion_reads_capture", "captured_bam", readset.sample.name, "captured.bam")
            out_dir = os.path.join("fusion_reads_capture", "realigned_bam", readset.sample.name)
            out_bam = os.path.join(out_dir, "realigned.bam")

            realign_job = extract_captured_reads_and_realn.extract_captured_reads_and_realn(captured_bam, out_bam, ini_section='extract_captured_reads_and_realn')

            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                realign_job
            ], name="extract_captured_reads_and_realn."+readset.sample.name)
            jobs.append(job)
        return jobs

    def filter_caputred_reads(self):
        """
        Compare capture alignment and realignment, print filtered fusion reads alignment and count summary                   
        """

        jobs = []
        for readset in self.readsets:
            sample_info_file = os.path.abspath(self.args.sampleinfo.name)                                                    

            out_dir = os.path.join("fusion_reads_capture", "filtered_result", readset.sample.name)
            captured_bam = os.path.join("fusion_reads_capture", "captured_bam", readset.sample.name, "captured.bam")
            realigned_bam = os.path.join("fusion_reads_capture", "realigned_bam", readset.sample.name, "realigned.bam")
            if readset.cram:
                filter_job = filter_caputred_reads.filter_caputred_reads(captured_bam, realigned_bam, os.path.basename(readset.cram), sample_info_file, out_dir, ini_section='filter_caputred_reads')
            elif readset.fastq1 or readset.bam:
                filter_job = filter_caputred_reads.filter_caputred_reads(captured_bam, realigned_bam, readset.sample.name, sample_info_file, out_dir, ini_section='filter_caputred_reads')
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                filter_job
            ], name="filter_caputred_reads."+readset.sample.name)
            jobs.append(job)
        return jobs

    def merge_and_summary(self):
        """
        Merge the results, extract summaries and get the statistics                                                          
        """

        jobs = []
        out_dir = os.path.join("fusion_reads_capture", "merged_summary_stats")
        sum_files = [os.path.join("fusion_reads_capture", "filtered_result", sample.name, "aln_and_summary") for sample in self.samples]
        sum_job = fusion_reads_capture_results_sum.merge_and_summary(sum_files, out_dir)
        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            sum_job
        ], name="merge_and_summary")

        jobs.append(job)
        return jobs

    def valfilter_cff_and_sample_enrichment(self):
        """
        Filter .cff file using valfilterX, where X is number of captured split reads, 
        and is specified as command line argument to pipeline. 
        Perform sample enrichment: Add cff entries for fusions found in additional samples 
        beyond the first sample in which the fusion was called 
        """
        jobs = []
        num_captured_reads = config.param('valfilter_cff_and_sample_enrichment', 'num_captured_reads', type='int') 
        merged_summary_dir = os.path.join("fusion_reads_capture", "merged_summary_stats")
        sum_file = os.path.join(merged_summary_dir , "merged.summary")
        #1) filter merged.summary file on num_captured_reads
        filter_merged_summary_job = valfilter.filter_merged_summary_file(sum_file, num_captured_reads)
        jobs.append(filter_merged_summary_job)
        #2) filter .cff and sample enrichment
        seq_len = config.param('repeat_filter', 'seq_len', type='int')
        #cff_file = os.path.join("fusions", "cff", "merged.cff.reann.dnasupp.bwafilter." + str(seq_len))
        cff_file = os.path.join("fusions", "cff", "merged.cff.reann.dnasupp")
        filtered_sum_file = os.path.join(merged_summary_dir, "merged." + str(num_captured_reads) +  ".summary")
        filter_cff_sample_enrichment_job = valfilter.filter_cff_and_sample_enrichment(filtered_sum_file, cff_file)
        jobs.append(filter_cff_sample_enrichment_job)
        return jobs

##END fusion reads capture pipeline

    def cluster_reann_dnasupp_file(self):
        """
        Reannotate DNA support (pair clusters) file. This step generates the final category/cluster file,
        merged.cff.reann.dnasupp.bwafilter.<seq-len>.cluster, which has the following columns:
            cluster_type, gene1, gene2, max_split_cnt, max_span_cnt, sample_type, disease, tools, inferred_fusion_type,
            gene1_on_bnd, gene1_close_to_bnd, gene2_on_bnd, gene2_close_to_bnd, dna_supp, samples
  
        Gene_Cluster    GTF2I   PGS1    2   -1  Tumor   LIS fusionmap   GeneFusion  True    True    True    True    -1  LIS_S8_L001

        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        cluster_job = merge_and_reannotate_cff_fusion.cluster_reann_dnasupp_file(out_dir)
        
        job = concat_jobs([
            cluster_job    
        ], name="cluster_reann_dnasupp_file")

        jobs.append(job)
        return jobs

    def custom_filters(self):
        """
        Includes blacklist_filter, callerfilter2 and ReadThrough filter 
        """
        jobs = []
        cff_dir = os.path.join("fusions", "cff")
        cluster_file=os.path.join(cff_dir, "merged.cff.renamed.reann.cluster")
        custom_filter_job= custom_filters.custom_filters(cff_dir, cluster_file)
        
        #category_filter_job = fusion_stats.category_filter(cluster_file)
        #caller_filter_job = fusion_stats.caller_filter(cluster_file)

        #covfilter_job = fusion_stats.covfilter(out_dir)
        filter_job = concat_jobs([
           custom_filter_job
        ], name="custom_filters")
        jobs.append(filter_job)
        return jobs
  
    def fusion_stats(self):
        """
        Outputs count files and plots about the detected gene fusions.
        """
        jobs = []
        cff_dir = os.path.join("fusions", "cff")
        out_dir= os.path.join("fusions", "fusion_stats")
        sampleinfo_file = os.path.relpath(self.args.sampleinfo.name, self.output_dir)

        fusion_stats_job = fusion_stats.fusion_stats(cff_dir, out_dir, sampleinfo_file)
        category_table_job = fusion_stats.generate_category_count_table(cff_dir, out_dir)
        category_barplot_job = fusion_stats.generate_categories_barplot(fusion_stats_dir=out_dir)
        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            fusion_stats_job,
            category_table_job,
            category_barplot_job
        ], name="fusion_stats")

        jobs.append(job)
        return jobs

    def validate_fusions(self):
        """
        Compares the pipeline output in merged.cff.reann.dnasupp.bwafilter.30.cluster with the predetermined
        fusion gene test file.This step should be run only with a test .bam/.fastq file, in order to confirm 
        detection of validated fusions which are known to bepresent in the sample.
        Requires --valfile flag, with corresponding file containing gene pairs
        """
        cff_dir = os.path.join("fusions", "cff")
        cluster_file = os.path.join(cff_dir, "merged.cff.reann.dnasupp.bwafilter.30.cluster")
        out_dir = os.path.join("fusions", "validate_fusions")
        jobs = []
        intersect_breakpoints_job = validate_fusions.intersect_breakpoints_with_bedfile(cluster_file, cff_dir)
        validate_fusions_job = validate_fusions.validate_fusions(cff_dir, out_dir, self.args.valfile.name)

        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            intersect_breakpoints_job,
            validate_fusions_job
        ], name="validate_fusions")

        jobs.append(job)
        return jobs

    def delete_fastqs(self):
        """
        Delete fastqs when all callers' jobs are finished
        """
        jobs = []
        for sample in self.samples:
            defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
            fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
            ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
            integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")
            result_file_list = [defuse_result, fusionmap_result, ericscript_result, integrate_result]
            del_job = delete_fastqs.delete_fastqs(sample.name, result_file_list)
            job = concat_jobs([
                Job(command="mkdir -p delete_fastqs"),
                del_job
            ], name="delete_fastqs." + sample.name)
            jobs.append(job)
        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.gunzip_fastq,
            self.merge_fastq,
            self.run_star_seqr,
            self.run_arriba,
            self.star_fusion,
            self.defuse,
            self.fusionmap,
            self.ericscript,
            self.tophat2,
            self.integrate,
            self.integrate_make_result_file,
            self.convert_fusion_results_to_cff,
            #self.fusioninspector
            #self.filter_cff_calls_using_fusioninspector_results,
            self.merge_cff_fusion,
            self.rename_genes,
            self.reannotate_cff_fusion,
            #self.check_dna_support_before_next_exon,
            #self.repeat_filter,
            # VALIDATION PIPELINE
            #self.build_fusion_and_head_gene_ref,
            #self.fastq_conversion_and_reads_capture,
            #self.extract_captured_reads_and_realn,
            #self.filter_caputred_reads,
            #self.merge_and_summary,
            #self.valfilter_cff_and_sample_enrichment,
            # MERGE .cff ENTRIES
            self.cluster_reann_dnasupp_file,
            self.custom_filters,
            #self.fusioninspector,
            #self.fusion_stats,
            #self.validate_fusions,
            self.delete_fastqs
        ]

if __name__ == '__main__':
    RnaFusion()
