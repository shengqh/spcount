import re
import pandas as pd

from .common_util import readFileMap

def krona(logger, treeFile, groupFile, outputPrefix):
  logger.info("Start krona ...")

  groupFileMap = readFileMap(groupFile)
  print(groupFileMap)

  tree_data = pd.read_csv(treeFile)
  print(tree_data.columns)

  logger.info("done")
