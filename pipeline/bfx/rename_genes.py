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

def rename_cff_file_genes(sample, tool, out_dir, ini_section='rename_genes'):

    cff_file = os.path.join(out_dir, sample.name + "." + tool + ".cff")
    cff_file_renamed = os.path.join(out_dir, sample.name + "." + tool + ".cff.renamed")

    return Job(
        [cff_file],
        [cff_file_renamed],
        [],
        command="""\
pwd; /hpf/largeprojects/ccmbio/mapostolides/mugqic_tools-my-version/python-tools/fusiontools/0.1.0/bin/validation_pipeline/Pipeline-scripts/rename_cff_file_genes-GENAP.py \\
  {cff_file} \\
  > {cff_file_renamed} """.format(
        cff_file=cff_file,
        cff_file_renamed=cff_file_renamed
        ),
        removable_files=[]
    )

