from context import spcount, taxonomy_file, refseq_file, output_dir

from spcount.database_util import prepare_segment_database

import logging

logger = logging.getLogger('assembly')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
prepare_segment_database(logger, taxonomy_file, refseq_file, 4751, output_dir, "20221124_fungi", 500, True)

