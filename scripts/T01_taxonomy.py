from context import spcount, taxonomy_file

from spcount.taxonomy_util import prepare_taxonomy

import logging

logger = logging.getLogger('taxonomy')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
prepare_taxonomy(logger, taxonomy_file)

