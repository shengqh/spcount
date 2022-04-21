from context import spcount, refseq_file, species_file, taxonomy_file

from spcount.count_util import bowtie_count

import logging

logger = logging.getLogger('assembly')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

input_list_file='/scratch/vickers_lab/projects/20220417_bacteria_genome/intermediate_data/refseq_bacteria_bowtie_count/result/WT_NSB_1.list'
count_file='/scratch/vickers_lab/projects/20200625_4893_2_RA_smRNA_mouse_v5_byTiger/intermediate_data/bowtie1_genome_unmapped_reads/result/WT_NSB_1_clipped_identical.unmapped.fastq.dupcount'
output_file="/scratch/vickers_lab/projects/20220417_bacteria_genome/intermediate_data/refseq_bacteria_bowtie_count/result/WT_NSB_1.txt.gz"

bowtie_count(logger, input_list_file, output_file, count_file, species_file, species_column='species')
