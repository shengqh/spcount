import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import spcount

taxonomy_file="/data/cqs/references/spcount/20220406_taxonomy.txt"
refseq_file="/data/cqs/references/spcount/20220406_assembly_summary_refseq.txt"
species_file="/data/cqs/references/spcount/20220406_bacteria.taxonomy.txt"