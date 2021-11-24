from context import spcount

from spcount.taxonomy_util import prepare_taxonomy

import logging

logger = logging.getLogger('taxonomy')
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
prepare_taxonomy(logger, "/data/cqs/references/spcount/20211118_taxonomy.txt")

