import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import spcount

#taxonomy_file="/data/cqs/references/spcount/20220406_taxonomy.txt"
#refseq_file="/data/cqs/references/spcount/20220406_assembly_summary_refseq.txt"
#species_file="/data/cqs/references/spcount/20220406_bacteria.taxonomy.txt"

taxonomy_file="/data1/shengq2/references/spcount/20221124.taxonomy.txt"
refseq_file="/data1/shengq2/references/spcount/20221124_assembly_summary.txt"
species_file="/data1/shengq2/references/spcount/20221124_virus.taxonomy.txt"
output_dir="/data1/shengq2/references/spcount"