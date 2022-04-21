from context import spcount, refseq_file, species_file, taxonomy_file

from spcount.count_util import count_table

import logging

logger = logging.getLogger('table')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')

input_list_file='/scratch/vickers_lab/projects/20220417_bacteria_genome/nonhost_genome/refseq_bacteria_table/result/RA_4893_2__fileList1.list'
output_prefix='/scratch/vickers_lab/projects/20220417_bacteria_genome/nonhost_genome/refseq_bacteria_table/result/RA_4893_2'

count_table(logger, input_list_file, output_prefix, species_file, species_column='species')
