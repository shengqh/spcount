from context import spcount, refseq_file

from spcount.database_util import download_assembly_summary

import logging

logger = logging.getLogger('assembly')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
download_assembly_summary(logger, refseq_file, "refseq")

