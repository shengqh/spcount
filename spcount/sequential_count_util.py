import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from _ctypes import ArgumentError

from .common_util import readFileMap
from .BowtieCountItem import BowtieCountItem, readBowtieTextFile

def prepare_bowtie_index(logger, sequence, thread):

  bowtieIndecies = readBowtieIndexList(categoryFile)

  if not slurm:
    for bowtieIndex in bowtieIndecies:
      logger.info("Building index for %s ..." % bowtieIndex.Fasta)
      indexDone = bowtieIndex.Index + ".index.done"
      if not os.path.exists(indexDone) or force:
        subprocess.call(['bowtie-build', '-q', '--threads', str(thread), bowtieIndex.Fasta, bowtieIndex.Index])
        open(indexDone, 'wt').close()
    return

def sequential_count(logger, fastqListFile, dbListFile, outputFile, countListFile):
  logger.info("Start count ...")

  fastqFileMap = readFileMap(fastqListFile)
  print(fastqFileMap)

  dbFileMap = readFileMap(dbListFile)
  print(dbFileMap)

  #print(bowtieFileMap)

  countFileMap = readFileMap(countListFile) if countListFile != None else {}
  print(countFileMap)
  
