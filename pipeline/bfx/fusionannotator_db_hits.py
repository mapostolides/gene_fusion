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

def fusionannotator_db_hits(outdir, cluster, ini_section='fusion_stats'):
    """ 
    Run FusionAnnotator and create cancer DB files
    """
    #output_file=    --> file which is "final output file" to check that step finished properly
    cancer_cluster = os.path.join(outdir, "merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter.CANCER_FUSIONS")
    normal_cluster = os.path.join(outdir, "merged.cff.renamed.reann.cluster.blck_filter.RT_filter.callerfilter2.ANC_filter..NORMALS")

    return Job(
        [cluster],
        [cancer_cluster, normal_cluster],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
module load libdb/4.7; source /home/mapostolides/miniconda3/etc/profile.d/conda.sh; /hpf/largeprojects/ccmbio/mapostolides/MODULES/RUN_BENCHMARKING_TOOLKIT/run_FusionAnnotator_cluster.sh {outdir} {cluster} """.format(
        outdir=outdir,
        cluster=cluster
        ),
        removable_files=[]
    )
