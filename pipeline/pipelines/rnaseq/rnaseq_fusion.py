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

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

from bfx import picard
from pipelines import common

from bfx import samtools_1_1
from bfx import defuse
from bfx import fusionmap
from bfx import tophat2
from bfx import integrate
from bfx import ericscript
from bfx import gunzip
from bfx import merge_fastq
from bfx import cff_conversion
from bfx import check_dna_support_before_next_exon
from bfx import merge_and_reannotate_cff_fusion
from bfx import repeat_filter
from bfx import delete_fastqs

import utils

log = logging.getLogger(__name__)

class RnaFusion(common.Illumina):
	"""
	RNAFusion Pipeline
	================
    The Gene Fusion pipeline identifies gene fusion events using RNA-seq FASTQ files.  
    
    Four separate tools detect fusion events: 
    [deFuse](https://sourceforge.net/p/defuse/wiki/DeFuse/), 
    [FusionMap](http://www.arrayserver.com/wiki/index.php?title=FusionMap), 
    [EricScript](https://sites.google.com/site/bioericscript/home), 
    and [INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home/).
    
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

    Example command: 
    
    python rnaseq_fusion.py -r readset.tsv -s 1-14 --sampleinfo disease.sampleinfo --dnabam disease.bam -c rnaseq.fusion.ini
    
	"""

	def __init__(self):
		# Add pipeline specific arguments
		self.argparser.add_argument("--sampleinfo", help="sample info file", type=file)
		self.argparser.add_argument("--dnabam", help="DNA bam list", type=file)
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
			trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
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


	def defuse(self):
		"""
		Run Defuse to call gene fusion
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
		Run FusionMap to call gene fusion
		"""
		jobs = []
		for sample in self.samples:
			fastq1, fastq2 = self.select_input_fastq(sample)
			out_dir = os.path.join("fusions", "fusionmap", sample.name)
			fusionmap_job = fusionmap.fusionmap(fastq1, fastq2, out_dir)
			job = concat_jobs([
				Job(command="mkdir -p " + out_dir),
				fusionmap_job,
				Job(command="ls " + out_dir + "/02_RNA*")
			
			], name="fusionmap." + sample.name)

			jobs.append(job)

		return jobs


	def ericscript(self):
		"""
		Run EricScript to call gene fusion
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
		Run Integrate to call gene fusion
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
			defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
			fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
			ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
			integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")

			tool_results = [("defuse", defuse_result), ("fusionmap", fusionmap_result), ("ericscript", ericscript_result), ("integrate", integrate_result)]
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
			for tool, result_file in tool_results:
				job = cff_conversion.cff_convert(sample.name, result_file, sampleinfo_file, tool, out_dir)
				job.command = job.command.strip()
				job_list.append(job)
		job = concat_jobs(job_list, name="cff_conversion")
		jobs.append(job)
		return jobs

	def merge_and_reannotate_cff_fusion(self):
		"""
		Merge all cff files into one single file and reannotate it with given annotation files
		"""
		jobs = []
		cff_files = []
		cff_dir = os.path.join("fusions", "cff")
		out_dir = os.path.join("fusions", "cff")
		tool_list = ["defuse", "fusionmap", "ericscript", "integrate"]
		for tool in tool_list:
			cff_files.extend([os.path.join(cff_dir, sample.name+"."+tool+".cff") for sample in self.samples])
		
		reann_job = merge_and_reannotate_cff_fusion.merge_and_reannotate_cff_fusion(cff_files, out_dir)
		
		job = concat_jobs([
			reann_job	
		], name="merge_and_reannotate_cff_fusion")

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

	def cluster_reann_dnasupp_file(self):
		"""
		Reannotate DNA support (pair clusters) file
		"""
		jobs = []
		out_dir = os.path.join("fusions", "cff")
		cluster_job = merge_and_reannotate_cff_fusion.cluster_reann_dnasupp_file(out_dir)
		
		job = concat_jobs([
			cluster_job	
		], name="cluster_reann_dnasupp_file")

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
			self.defuse,
			self.fusionmap,
			self.ericscript,
			self.tophat2,
			self.integrate,
			self.integrate_make_result_file,
			self.convert_fusion_results_to_cff,
			self.merge_and_reannotate_cff_fusion,
			self.check_dna_support_before_next_exon,
			self.repeat_filter,
			self.cluster_reann_dnasupp_file,
			self.delete_fastqs
		]

if __name__ == '__main__':
	RnaFusion()
