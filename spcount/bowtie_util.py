import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from _ctypes import ArgumentError

from .BowtieIndex import BowtieIndexItem, readBowtieIndexList

def bowtie(logger, inputFile, outputFile, databaseListFile, thread, isFasta=False):
  logger.info("Start bowtie ...")
  bowtieIndexList = readBowtieIndexList(databaseListFile)

  for bowtieIndex in bowtieIndexList:
    if not os.path.exists(bowtieIndex.Index + ".rev.2.ebwt") and not os.path.exists(bowtieIndex.Index + ".rev.2.ebwtl") :
      raise ArgumentError("Bowtie index not exists: %s" % bowtieIndex.Index)

  categories = sorted(set([bi.Category for bi in bowtieIndexList])) 
  logfile = outputFile + ".log"
  tmpFile = outputFile + ".tmp"
  with open(tmpFile, "wt") as fout:
    with open(logfile, "wt") as flog:
      for category in categories:
        logger.info("Searching category %s ..." % category)
        flog.write(">" + category + "\n")
        flog.flush()
        bowtieIndecies = [bi.Index for bi in bowtieIndexList if bi.Category == category]

        bowtieCount = 0

        inputFastq = inputFile
        for bowtieIndex in bowtieIndecies:
          bowtieCount = bowtieCount + 1
          logger.info("  Searching to %s ..." % bowtieIndex)
          outfile = "%s.%d.out" % (outputFile, bowtieCount)

          if isFasta:
            unmapped = "%s.%d.unmapped.fasta" % (outputFile, bowtieCount)
            subprocess.call(['bowtie', '-f', '-k', '1', '-v', '0', '-p', str(thread), '--no-unal', '--un', unmapped, bowtieIndex, inputFastq, outfile], 
              stderr=flog)
          else:
            unmapped = "%s.%d.unmapped.fastq" % (outputFile, bowtieCount)
            subprocess.call(['bowtie', '-k', '1', '-v', '0', '-p', str(thread), '--no-unal', '--un', unmapped, bowtieIndex, inputFastq, outfile], 
              stderr=flog)

          with open(outfile, "rt") as fin:
            for line in fin:
              parts = line.rstrip().split('\t')
              fout.write("%s\t%s\t%s\t%s\t%s\n" % (parts[0].split(' ')[0], parts[1], parts[2], parts[3], category))

          os.remove(outfile)
          
          if bowtieCount > 1:
            os.remove(inputFastq)
          inputFastq = unmapped
        
        if inputFastq != inputFile:
          os.remove(inputFastq)
  
  os.rename(tmpFile, outputFile)

  logger.info("done")

def bowtie_fastq2fasta(logger, inputFile, outputFile, databaseListFile, thread):
  logger.info("Start bowtie_fastq2fasta ...")
  bowtieIndexList = readBowtieIndexList(databaseListFile)

  for bowtieIndex in bowtieIndexList:
    if not os.path.exists(bowtieIndex.Index + ".rev.2.ebwt") and not os.path.exists(bowtieIndex.Index + ".rev.2.ebwtl") :
      raise ArgumentError("Bowtie index not exists: %s" % bowtieIndex.Index)

  inputFasta = outputFile + ".fasta"
  fin = gzip.open(inputFile, "rt") if inputFile.endswith(".gz") else open(inputFile, "rt")
  with fin:
    with open(inputFasta, "wt") as fout:
      while True:
        query = fin.readline()
        if not query:
          break
        sequence = fin.readline().rstrip()
        fin.readline()
        fin.readline()
        queryName = query.rstrip().split('\t')[0].split(' ')[0]
        fout.write(">%s\n%s\n" %(queryName[1:], sequence))
  
  bowtie(logger, inputFasta, outputFile, databaseListFile, thread, True)
  
if __name__ == "__main__":
  logger = logging.getLogger('sequenceBowtie')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  bowtie(logger, 
            "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.gz", 
            "/scratch/cqs/shengq2/temp/VLDL_WZ_bacteria.txt",
            "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.slim.list",
            32)
