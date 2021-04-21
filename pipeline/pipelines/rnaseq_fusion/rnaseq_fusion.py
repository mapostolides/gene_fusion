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
#callers
from bfx import star_seqr 
from bfx import arriba
from bfx import star_fusion
from bfx import defuse
from bfx import fusionmap
from bfx import tophat2
from bfx import integrate
from bfx import ericscript
from bfx import chimerascan 
from bfx import metafusion 
from bfx import metafusion_isohunter 
from bfx import metafusion_clinical 

from bfx import gunzip
from bfx import merge_fastq
from bfx import cff_conversion
from bfx import merge_and_reannotate_cff_fusion
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
        #self.argparser.add_argument("--dnabam", help="DNA bam list", type=file)
        # add optional fusion validation file for pipeline validation mode
        self.argparser.add_argument("--valfile", required=False, help="fusion validation set file", type=file)
        self.argparser.add_argument("--callers", required=False, help="comma-separated string of callers to be used")
        self.argparser.add_argument("--keep_bams", required=False,  action='store_true', help="if this flag is used, keep caller bams")
        self.argparser.add_argument("--database", required=False, help="specifies database path for historical fusions, FP list and clinical fusions")
        #Class variables
        self.tool_list = ["star_seqr", "arriba", "star_fusion", "fusionmap", "ericscript", "integrate", "defuse"]
        #self.tool_list = ["star_seqr", "arriba", "star_fusion", "fusionmap", "ericscript", "defuse"]
        #self.tool_list = ["arriba", "star_fusion", "fusionmap", "ericscript", "integrate", "defuse"]
        super(RnaFusion, self).__init__()


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

    def chimerascan(self):
        """
        Run chimerascan to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "chimerascan", sample.name)
            chimerascan_job = chimerascan.run(fastq1, fastq2, out_dir )
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                Job(command="rm -r " + out_dir),
                chimerascan_job
            ], name="chimerascan." + sample.name)

            jobs.append(job)

        return jobs

    def run_arriba(self):
        """
        """

        jobs = []
        for sample in self.samples:
            if len(sample.readsets) > 1:
                raise Exception("Error: only one read set per sample allowed")
            if sample.readsets[0].bam:#.bam input
                fastq_dir = os.path.join("fusions", "picard_sam_to_fastq", sample.name)
                bam = sample.readsets[0].bam
                left_fastq=os.path.join(self._output_dir, fastq_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", bam)))
                right_fastq=os.path.join(self._output_dir, fastq_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", bam)))
            elif sample.readsets[0].fastq2 and sample.readsets[0].fastq2.split(".")[-1] == "gz":
                left_fastq = sample.readsets[0].fastq1
                right_fastq = sample.readsets[0].fastq2
            else:
                raise Exception("Error: only .bam and .fastq.gz inputs allowed")
            output_dir = os.path.join("fusions", "arriba", sample.name)
            #JOBS
            chgdir_job = Job(command="cd " + output_dir)
            back_to_outdir_job = Job(command="cd " + self._output_dir)
            #CONCAT
            job = concat_jobs([
                Job(command="mkdir -p " + output_dir),
                chgdir_job,
                arriba.run(left_fastq, right_fastq, self._output_dir, output_dir, keep_bam=self.args.keep_bams),
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
                bam = sample.readsets[0].bam
                #fastq1 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", out_bam)))
                #fastq2 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", out_bam)))
                left_fastq=os.path.join(fastq_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", bam)))
                right_fastq=os.path.join(fastq_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", bam)))
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
                star_seqr.run(left_fastq, right_fastq, output_dir, sample.name, keep_bam=self.args.keep_bams)
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
            star_fusion_job = star_fusion.star_fusion(fastq1, fastq2, CTAT_resource_lib, out_dir, keep_bam=self.args.keep_bams)
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
            defuse_job = defuse.defuse(fastq1, fastq2, out_dir, keep_bam=self.args.keep_bams)
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
            ericscript_job = ericscript.ericscript(fastq1, fastq2, out_dir, keep_bam=self.args.keep_bams)
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
            
            # Define result files
            #output_file = os.path.join(output_dir, prefix + "_STAR-SEQR", prefix  + "_STAR-SEQR_candidates.txt")
            #star_seqr_result = os.path.join("fusions", "star_seqr", sample.name, "out_STAR-SEQR", "out_STAR-SEQR_candidates.txt")
            star_seqr_result = os.path.join("fusions", "star_seqr", sample.name, "out_STAR-SEQR_candidates.txt")
            #print >> sys.stderr, star_seqr_result 
            arriba_result = os.path.join("fusions", "arriba", sample.name, "fusions.tsv")
            #star_fusion_result = os.path.join("fusions", "star_fusion", sample.name, "star-fusion.fusion_predictions.abridged.tsv")
            star_fusion_result = os.path.join("fusions", "star_fusion", sample.name, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
            defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
            fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
            ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
            integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")
            # Build tool_results list based on self.tool_list
            result_file_dict = {"star_seqr": star_seqr_result, "arriba": arriba_result, 
                                "star_fusion": star_fusion_result, "defuse": defuse_result, 
                                "fusionmap": fusionmap_result, "ericscript": ericscript_result, 
                                "integrate": integrate_result
                                }
            tool_results = [(key, result_file_dict[key]) for key in result_file_dict.keys() if key in self.tool_list]
            #tool_results = [("star_seqr",star_seqr_result), ("arriba", arriba_result), ("star_fusion", star_fusion_result), ("defuse", defuse_result), ("fusionmap", fusionmap_result), ("ericscript", ericscript_result), ("integrate", integrate_result)]
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
            cff_files.extend([os.path.join(cff_dir, sample.name + "." + tool + ".cff") for sample in self.samples])
        merge_job = merge_and_reannotate_cff_fusion.merge_cff_fusion(cff_files, out_dir)
        
        job = concat_jobs([ merge_job ], name="merge_cff_fusion")
        jobs.append(job)
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

    def MetaFusion(self):
        """
        Run MetaFusion
        """
        jobs = []
        cff_dir_abspath = os.path.join(self._output_dir, "fusions", "cff")
        out_dir_abspath = os.path.join(self._output_dir, "fusions", "metafusion")
        metafusion_job = metafusion.run_metafusion_singularity(out_dir_abspath)
        #metafusion_job.name = "MetaFusion"
        job = concat_jobs([
            Job(command="mkdir -p " + out_dir_abspath),
            metafusion_job 
        ], name="MetaFusion")

        jobs.append(job)
        
        return jobs


    def MetaFusion_IsoHunter(self):
        """
        Run MetaFusion.IsoHunter
        """
        jobs = []
        out_dir_abspath = self._output_dir
        isohunter_outdir = os.path.join("fusions", "metafusion_isohunter") 
        metafusion_job = metafusion_isohunter.run_isohunter_singularity(out_dir_abspath)
        job = concat_jobs([
            Job(command="mkdir -p " + isohunter_outdir),
            metafusion_job
        ], name="MetaFusion.IsoHunter")

        jobs.append(job)

        return jobs

    def MetaFusion_clinical(self):
        """
        Run MetaFusion.IsoHunter.clinical
        """
        jobs = []
        out_dir_abspath = self._output_dir
        metafusion_outdir = os.path.join("fusions", "metafusion_clinical")
        metafusion_job = metafusion_clinical.run_metafusion_clinical(out_dir_abspath, self.args.database)
        job = concat_jobs([
            Job(command="mkdir -p " + metafusion_outdir),
            metafusion_job
        ], name="MetaFusion.clinical")

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
            star_seqr_result = os.path.join("fusions", "star_seqr", sample.name, "out_STAR-SEQR_candidates.txt")
            arriba_result = os.path.join("fusions", "arriba", sample.name, "fusions.tsv")
            star_fusion_result = os.path.join("fusions", "star_fusion", sample.name, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")

            #result_file_list = [defuse_result, fusionmap_result, ericscript_result, integrate_result, star_seqr_result, arriba_result, star_fusion_result] 
            result_file_list = [defuse_result, fusionmap_result]
            del_job = delete_fastqs.delete_fastqs(sample.name, result_file_list) 
            job = concat_jobs([
                Job(command="mkdir -p delete_fastqs"),
                del_job
            ], name="delete_fastqs." + sample.name)
            #job = concat_jobs([
            #    Job(command="mkdir -p delete_fastqs")
            #], name="delete_fastqs." + sample.name)
            jobs.append(job)
            # DELETE BAMS JOB (one across all samples)
        del_bams_job = concat_jobs([
               delete_fastqs.delete_bams(result_file_list, self._output_dir)
               ], name="delete_bams") 
        jobs.append(del_bams_job)
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
            self.merge_cff_fusion,
    #        self.MetaFusion,
            self.MetaFusion_clinical,
            self.delete_fastqs,
            self.MetaFusion_IsoHunter,
            self.chimerascan
        ]

if __name__ == '__main__':
    RnaFusion()
