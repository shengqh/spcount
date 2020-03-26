import argparse
import sys
import logging
import os
import csv
import ftplib
import gzip
import subprocess
from datetime import datetime
from io import StringIO
from .BowtieIndex import BowtieIndexItem, readBowtieIndexList, writeBowtieIndexList
from .Taxonomy import TaxonomyItem, TaxonomyTree
from .Category import CategoryItem

def getCategory(source):
  return(source.replace(" group", "").replace("/", "_").replace(" Bacteria", "").replace(" ", "_"))

def combine_category_fasta_file(logger, maxGenomeInFile, localDir, prefix, targetFile, localFileMap):
  categoryFile = targetFile + ".list"
  categoryFileDone = categoryFile + ".done"

  if os.path.exists(categoryFileDone):
    return(categoryFile)

  totalCount = len(localFileMap)

  bowtieIndecies = []
  foutMap = {}
  currentCount = 0
  with open(targetFile, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')

      category = getCategory(parts[2])

      if not category in foutMap.keys():
        subName = "%s%s_%03d" % (prefix, category, 1)
        bowtieIndex = os.path.join(localDir, subName)
        categoryFastaFile = bowtieIndex + ".fasta"
        logger.info("Init " + categoryFastaFile)
        #cat = CategoryItem(None, 1, 1)
        cat = CategoryItem(open(categoryFastaFile, "wb"), 1, 1)
        foutMap[category] = cat
        bowtieIndecies.append(BowtieIndexItem(bowtieIndex, category, categoryFastaFile))
      else:
        cat = foutMap[category]
        if cat.NumberOfGenome >= maxGenomeInFile:
          cat.FileHandle.close()
          cat.Index = cat.Index + 1
          cat.NumberOfGenome = 1
          subName = "%s%s_%03d" % (prefix, category, cat.Index)
          bowtieIndex = os.path.join(localDir, subName)
          categoryFastaFile = bowtieIndex + ".fasta"
          logger.info("Re-init " + categoryFastaFile)
          cat.FileHandle = open(categoryFastaFile, "wb")
          bowtieIndecies.append(BowtieIndexItem(bowtieIndex, category, categoryFastaFile))
        else:
          cat.NumberOfGenome = cat.NumberOfGenome + 1

      categoryOut = cat.FileHandle
      
      localFnaFile = localFileMap[parts[3]]
      
      currentCount = currentCount + 1
      logger.info("Merge %d/%d: %s ..." % (currentCount, totalCount, localFnaFile))
      with gzip.open(localFnaFile, 'rb') as f:
        file_content = f.read()
        categoryOut.write(file_content)
      
  for fout in foutMap.values():
    fout.FileHandle.close()

  bowtieIndecies.sort(key=lambda x: x.Fasta)
  writeBowtieIndexList(categoryFile, bowtieIndecies)

  open(categoryFileDone, 'wt').close()
  return(categoryFile)

def open_ftp(name):
  rootDir = "/genomes/refseq/" + name.lower() + "/"
  result = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
  result.login("anonymous", "quanhu.sheng.1@vumc.org")
  result.cwd(rootDir)
  return(result)

def download_assembly_summary(logger, name, localDir, date_time):
  summaryFile = "assembly_summary.txt"
  result = os.path.join(localDir, date_time + summaryFile)
  
  if os.path.exists(result):
    return(result)

  logger.info("Downloading assembly_summary.txt ...")
  with open_ftp(name) as ftp:
    with open(result, "wb") as f:
      ftp.retrbinary("RETR " + summaryFile, f.write)

  return(result)

def extract_complete_genome(logger, rootFile, idMap):
  result = rootFile + ".files"

  if os.path.exists(result):
    return(result)

  with open(result, "wt") as fout:
    fout.write("Type\tAccession\tCategory\tFile\n")
    logger.info("Reading assembly_summary.txt ...")
    with open(rootFile, "rt") as fin:
      count = 0
      for line in fin:
        if line.startswith('#'):
          continue

        parts = line.split('\t')
        taxid = parts[5]
  
        if taxid not in idMap:
          continue
  
        version_status = parts[10]
        assembly_level = parts[11]
  
        if (assembly_level != "Complete Genome") or (version_status != "latest"):
          continue
  
        count = count + 1
        if count % 1000 == 0:
          #break
          print(count)
  
        accession = parts[0]
        url = parts[19]
        folderName = os.path.basename(url)
        fnaFile = url + "/" + folderName + "_genomic.fna.gz"
        urlFile = fnaFile.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
        fout.write("Complete Genome\t%s\t%s\t%s\n" % (accession, idMap[taxid], urlFile))

  return(result)

def download_assembly_genomes(logger, name, ftpname, genomeRootDir, targetFile):
  logger.info("Checking %s ... " % targetFile)

  totalCount = 0
  with open(targetFile, "rt") as fin:
    fin.readline()
    for line in fin:
      totalCount = totalCount + 1

  logger.info("Downloading %d genomes to cache folder %s ..." % (totalCount, genomeRootDir))
  currentCount = 0
  localFileMap = {}

  ftp = None
  with open(targetFile, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')

      category = getCategory(parts[2])
      fnaFile = parts[3]

      categoryDir = os.path.join(genomeRootDir, category)

      if not os.path.exists(categoryDir):
        os.mkdir(categoryDir)

      currentCount = currentCount + 1

      tmpFile = os.path.join(categoryDir, os.path.basename(fnaFile))

      localFileMap[fnaFile] = tmpFile

      tmpDoneFile = tmpFile + ".done"
      if os.path.exists(tmpDoneFile):
        continue

      logger.info("Downloading %d/%d: %s to %s ..." % (currentCount, totalCount, fnaFile, tmpFile))
      with open(tmpFile, "wb") as f:
        if ftp == None:
          ftp = open_ftp(ftpname)
        ftp.retrbinary("RETR " + fnaFile, f.write)
      open(tmpDoneFile, 'wt').close()

  if (ftp != None):
    ftp.close()

  return(localFileMap)

def prepare_database(logger, taxonomyRootId, outputFolder, taxonomyFile, maxGenomeInFile, prefix):
  tree = TaxonomyTree()

  nameToFolder = {'Viruses':'viral'}

  logger.info("Reading taxonomy from %s ..." % taxonomyFile)
  tree.ReadFromFile(taxonomyFile)

  idMap = tree.BuildChildIdParentNameMap(taxonomyRootId)
  name = tree.GetItem(taxonomyRootId).Name

  if name in nameToFolder:
    ftpname = nameToFolder[name]
  else:
    ftpname = name

  logger.info("ftpname=%s" % ftpname)

  localDir = os.path.join(outputFolder, name.lower())
  if not os.path.exists(localDir):
    os.mkdir(localDir)

  cacheDir = os.path.join(localDir, "cache")
  if not os.path.exists(cacheDir):
    os.mkdir(cacheDir)

  fastaDir = os.path.join(localDir, "fasta")
  if not os.path.exists(fastaDir):
    os.mkdir(fastaDir)

  rootFile = download_assembly_summary(logger, ftpname, localDir, prefix)

  targetFile = extract_complete_genome(logger, rootFile, idMap)
  
  localFileMap = download_assembly_genomes(logger, name, ftpname, cacheDir, targetFile)

  categoryFile = combine_category_fasta_file(logger, maxGenomeInFile, fastaDir, prefix, targetFile, localFileMap)

  logger.info("Done.")

def prepare_index(logger, categoryFile, thread, force=False, slurm=False, slurmEmail=None, slurmTemplate=None):
  bowtieIndecies = readBowtieIndexList(categoryFile)

  if not slurm:
    for bowtieIndex in bowtieIndecies:
      logger.info("Building index for %s ..." % bowtieIndex.Fasta)
      indexDone = bowtieIndex.Index + ".index.done"
      if not os.path.exists(indexDone) or force:
        subprocess.call(['bowtie-build', '-q', '--threads', str(thread), bowtieIndex.Fasta, bowtieIndex.Index])
        open(indexDone, 'wt').close()
    return

  if slurmTemplate == None:
    slurmTemplate = os.path.join(os.path.dirname(os.path.realpath(__file__)), "slurm.template")

  if not os.path.exists(slurmTemplate):
    raise ArgumentError("Slurm template not exists: %s" % slurmTemplate)

  with open(slurmTemplate, "rt") as fin:
    slurmLines = [line.rstrip() for line in fin]

  slurmFolder = os.path.join(os.path.dirname(categoryFile), "slurm")
  if not os.path.exists(slurmFolder):
    os.mkdir(slurmFolder)

  submitFile = os.path.join(slurmFolder, "submit.sh")
  with open(submitFile, "wt") as fout:
    for bowtieIndex in bowtieIndecies:
      indexDone = bowtieIndex.Index + ".index.done"
      command = "bowtie-build -q --threads %s %s %s\necho '' > %s\n" % (thread, bowtieIndex.Fasta, bowtieIndex.Index, indexDone )

      slurmFile = os.path.join(slurmFolder, os.path.basename(bowtieIndex.Index + ".slurm"))
      with open(slurmFile, "wt") as fslurm:
        for line in slurmLines:
          if "__EMAIL__" in line:
            fslurm.write(line.replace("__EMAIL__", slurmEmail) + "\n")
          elif "__THREAD__" in line:
            fslurm.write(line.replace("__THREAD__", str(thread)) + "\n")
          elif "__COMMAND__" in line:
            fslurm.write(line.replace("__COMMAND__", command) + "\n")
          elif "__LOG__" in line:
            fslurm.write(line.replace("__LOG__", slurmFile + ".log") + "\n")
          else:
            fslurm.write(line + "\n")

      if not force:
        fout.write("if [[ ! -e %s ]]; then\n  sbatch %s\nfi\n\n" % (indexDone, slurmFile))
      else:
        fout.write("sbatch %s\n\n" % slurmFile )

  logger.info("Please submit job using %s" % submitFile)
  return
