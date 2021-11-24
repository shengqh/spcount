import argparse
import sys
import logging
import os
import csv
import ftplib
import gzip
import requests
import shutil
import tarfile
import subprocess
import time
from datetime import datetime
from io import BytesIO, StringIO, TextIOWrapper

from .GenomeItem import GenomeItem, writeGenomeItems, readGenomeItems
from .Taxonomy import TaxonomyItem, TaxonomyTree

class CategoryItem(object):
  def __init__(self, fileHandle, index, numberOfGenome):
    self.FileHandle = fileHandle
    self.Index = index
    self.NumberOfGenome = numberOfGenome

def extract_unique_genome_priority(logger, resultFile, databaseFile, taxonomyFile, taxonomyRootId, priorityLevels):
  logger.info("Reading taxonomy from %s ..." % taxonomyFile)
  tree = TaxonomyTree()
  tree.ReadFromFile(taxonomyFile)

  idCategoryMap = tree.BuildChildIdParentNameMap(taxonomyRootId)

  levels = priorityLevels[::-1]

  idItemMap = {}
  logger.info("Reading %s ..." % databaseFile)
  with open(databaseFile, "rt") as fin:
    for line in fin:
      if line.startswith('#'):
        continue

      parts = line.split('\t')
      version_status = parts[10]

      if (version_status != "latest"):
        continue

      url = parts[19]
      if url == 'na':
        continue

      taxidstr = parts[5]
      if taxidstr not in idCategoryMap:
        continue

      taxonomyId = int(taxidstr)

      assemblyLevel = parts[11]
      assemblyLevelIndex = levels.index(assemblyLevel)

      if taxonomyId in idItemMap:
        if idItemMap[taxonomyId].AssemblyLevelIndex > assemblyLevelIndex:
          continue

      accession = parts[0]
      name = parts[7]
        
      folderName = os.path.basename(url)
      fnaFile = url + "/" + folderName + "_genomic.fna.gz"
      urlFile = fnaFile.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
      category = getCategory(idCategoryMap[taxidstr])
      
      idItemMap[taxonomyId] = GenomeItem(category, taxonomyId, accession, name, assemblyLevel, urlFile, assemblyLevelIndex)

  giItems = [idItemMap[taxid] for taxid in sorted(idItemMap.keys())]
  writeGenomeItems(resultFile, giItems)

def extract_unique_genome(logger, databaseFile, taxonomyFile, taxonomyRootId, completeGenomeOnly=False):
  if completeGenomeOnly:
    priorityLevels = ["Complete Genome"]
  else:
    priorityLevels = ["Complete Genome", "Scaffold", "Contig", "Chromosome"]

  return(extract_unique_genome_priority(logger, databaseFile, taxonomyFile, taxonomyRootId, priorityLevels))
