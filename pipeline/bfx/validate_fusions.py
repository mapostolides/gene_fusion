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

def validate_fusions(outdir, cluster, truth_fusions, ini_section='fusion_stats'):
    """ 
    Intersects the breakpoints of called fusions in the file merged.cff.reann.dnasupp.bwafilter.30.cluster
    """
    cff=os.path.join(outdir, "merged.cff.renamed.reann")
    #output_file=    --> file which is "final output file" to check that step finished properly

    return Job(
        [cluster],
        [],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
module load libdb/4.7; /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/benchmarking_cluster-GENAP.sh {outdir} {truth_fusions} {cff} {cluster} true true true true false""".format(
        outdir=outdir,
        truth_fusions=truth_fusions,
        cff=cff,
        cluster=cluster
        ),
        removable_files=[]
    )
