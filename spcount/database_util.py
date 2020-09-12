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
from datetime import datetime
from io import BytesIO, StringIO, TextIOWrapper

from .BowtieIndex import BowtieIndexItem, readBowtieIndexList, writeBowtieIndexList
from .GenomeItem import GenomeItem, writeGenomeItems, readGenomeItems
from .Taxonomy import TaxonomyItem, TaxonomyTree

class CategoryItem(object):
  def __init__(self, fileHandle, index, numberOfGenome):
    self.FileHandle = fileHandle
    self.Index = index
    self.NumberOfGenome = numberOfGenome

def getCategory(source):
  return(source.replace(" group", "").replace("/", "_").replace(" Bacteria", "").replace(" ", "_"))

def prepare_taxonomy(logger, outputFile):
  if os.path.exists(outputFile):
    return

  logger.info("Preparing taxonomy file %s ..." % outputFile)

  #download taxonomy from ncbi
  url="https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  targetFile = outputFile + "." + os.path.basename(url)
  logger.info("Downloading %s ..." % url)
  r = requests.get(url, verify=False, stream=True)
  r.raw.decode_content = True
  with open(targetFile, 'wb') as f:
    shutil.copyfileobj(r.raw, f)

  #prepare taxonomy file with three columns: id, parentid, name
  with tarfile.open(targetFile, "r:gz") as tar:
    taxonomies = []
    logger.info("Reading names.dmp ...")
    with TextIOWrapper(tar.extractfile(tar.getmember('names.dmp'))) as fin:
      for line in fin:
        if "scientific name" in line:
          parts = line.split('|')
          id = parts[0].strip()
          name = parts[1].strip()
          taxonomies.append(TaxonomyItem(name, id))

    taxonomyMap = {item.Id:item for item in taxonomies}

    logger.info("Reading nodes.dmp ...")
    with TextIOWrapper(tar.extractfile(tar.getmember('nodes.dmp'))) as fin:
      for line in fin:
        parts = line.split('|')
        childId = parts[0].strip()
        parentId = parts[1].strip()
        taxonomyMap[childId].ParentId = parentId

    with open(outputFile, "wt") as fout:
      fout.write("Id\tParentId\tScientificName\n")
      for item in taxonomies:
        fout.write("%s\t%s\t%s\n" % (item.Id, item.ParentId, item.Name))
  
  os.remove(targetFile)
  return

def combine_category_fasta_file(logger, maxGenomeInFile, localDir, prefix, targetFile, localFileMap):
  categoryFile = targetFile + ".list"
  categoryFileDone = get_done_file(categoryFile)

  if os.path.exists(categoryFileDone):
    return(categoryFile)

  totalCount = len(localFileMap)

  bowtieIndecies = []
  foutMap = {}
  currentCount = 0

  giItems = readGenomeItems(targetFile)

  for gi in giItems:
    category = gi.Category

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
    
    localFnaFile = localFileMap[gi.UrlFile]
    
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

def combine_fasta_file(logger, localDir, prefix, targetFile, localFileMap):
  fasta_file = os.path.join(localDir, prefix + "genomes.fasta")
  map_file =  os.path.join(localDir, prefix + "genomes.map")

  giItems = readGenomeItems(targetFile)
  currentCount = 0
  totalCount = len(giItems)
  with open(fasta_file, "wb") as fFasta:
    with open(map_file, "wt") as fMap:
      fMap.write("Chromosome\tSpecies\tCategory\n")
      for gi in giItems:
        category = gi.Category
        name = gi.Name
        localFnaFile = localFileMap[gi.UrlFile]

        currentCount = currentCount + 1
        logger.info("Merge %d/%d: %s ..." % (currentCount, totalCount, localFnaFile))
        with gzip.open(localFnaFile, 'rb') as f:
          file_content = f.read()
          fFasta.write(file_content)
          with TextIOWrapper(BytesIO(file_content)) as fin:
            for line in fin:
              if line.startswith(">"):
                logger.info(line.rstrip())
                parts = line.rstrip().split(' ')
                chrom = parts[0]
                fMap.write("%s\t%s\t%s\n" % (chrom, name, category))
          
  return(fasta_file)

def open_ftp():
  rootDir = ""
  result = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
  result.login("anonymous", "quanhu.sheng.1@vumc.org")
  result.cwd(rootDir)
  return(result)

def download_assembly_summary(logger, localDir, prefix, database="genbank"):
  summaryFile =  "/genomes/%s/assembly_summary_%s.txt" % (database, database)
  result = os.path.join(localDir, prefix + os.path.basename(summaryFile))
  
  if os.path.exists(result):
    return(result)

  logger.info("Downloading %s ..." % summaryFile)
  with open_ftp() as ftp:
    with open(result, "wb") as f:
      ftp.retrbinary("RETR " + summaryFile, f.write)

  return(result)

def extract_unique_genome_priority(logger, rootFile, idCategoryMap, taxonomyRootId, priorityLevels):
  priorStr = "_".join([s[0:4] for s in priorityLevels])

  result = "%s.%s.%s.files" % (rootFile, taxonomyRootId, priorStr)

  levels = priorityLevels[::-1]

  idItemMap = {}
  logger.info("Reading %s ..." % rootFile)
  with open(rootFile, "rt") as fin:
    for line in fin:
      if line.startswith('#'):
        continue

      parts = line.split('\t')
      version_status = parts[10]

      if (version_status != "latest"):
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
      url = parts[19]

      if url == 'na':
        continue
        
      folderName = os.path.basename(url)
      fnaFile = url + "/" + folderName + "_genomic.fna.gz"
      urlFile = fnaFile.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
      category = getCategory(idCategoryMap[taxidstr])
      
      idItemMap[taxonomyId] = GenomeItem(category, taxonomyId, accession, name, assemblyLevel, urlFile, assemblyLevelIndex)

  giItems = [idItemMap[taxid] for taxid in sorted(idItemMap.keys())]
  writeGenomeItems(result, giItems)
  return(result)

def extract_unique_genome(logger, rootFile, idMap, taxonomyRootId, completeGenomeOnly=False):
  if completeGenomeOnly:
    priorityLevels = ["Complete Genome"]
  else:
    priorityLevels = ["Complete Genome", "Scaffold", "Contig", "Chromosome"]

  return(extract_unique_genome_priority(logger, rootFile, idMap, taxonomyRootId, priorityLevels))

def get_done_file(fileName):
  return(fileName + ".done")

def download_assembly_genomes(logger, genomeRootDir, targetFile):
  logger.info("Checking %s ... " % targetFile)

  giItems = readGenomeItems(targetFile)
  totalCount = len(giItems)
  logger.info("Total %d genomes to cache folder %s ..." % (totalCount, genomeRootDir))
  localFileMap = {}

  for gi in giItems:
    categoryDir = os.path.join(genomeRootDir, gi.Category)

    if not os.path.exists(categoryDir):
      os.mkdir(categoryDir)

    localFile = os.path.join(categoryDir, os.path.basename(gi.UrlFile))
    localFileMap[gi.UrlFile] = localFile

  waitingFileMap = {fnaFile:localFileMap[fnaFile] for fnaFile in localFileMap if not os.path.exists(get_done_file(localFileMap[fnaFile]))}
  totalCount = len(waitingFileMap)
  logger.info("Downloading %d genomes to cache folder %s ..." % (totalCount, genomeRootDir))

  if totalCount > 0:
    with open_ftp() as ftp:
      currentCount = 0
      for fnaFile in waitingFileMap.keys():
        currentCount = currentCount + 1
        localFile = waitingFileMap[fnaFile]
        localDoneFile = get_done_file(localFile)

        logger.info("Downloading %d/%d: %s to %s ..." % (currentCount, totalCount, fnaFile, localFile))
        with open(localFile, "wb") as f:
          ftp.retrbinary("RETR " + fnaFile, f.write)
        open(localDoneFile, 'wt').close()

  return(localFileMap)

def prepare_database(logger, taxonomyRootId, outputFolder, maxGenomeInFile, prefix, database="genbank"):
  taxonomyFile = os.path.join(outputFolder, prefix + "taxonomy.txt")
  
  prepare_taxonomy(logger, taxonomyFile)

  logger.info("Reading taxonomy from %s ..." % taxonomyFile)
  tree = TaxonomyTree()
  tree.ReadFromFile(taxonomyFile)

  idMap = tree.BuildChildIdParentNameMap(taxonomyRootId)
  name = tree.GetItem(taxonomyRootId).Name.replace(" ", "_") + "." + taxonomyRootId

  localDir = os.path.join(outputFolder, name.lower())
  if not os.path.exists(localDir):
    os.mkdir(localDir)

  cacheDir = os.path.join(localDir, "cache")
  if not os.path.exists(cacheDir):
    os.mkdir(cacheDir)

  fastaDir = os.path.join(localDir, "fasta")
  if not os.path.exists(fastaDir):
    os.mkdir(fastaDir)

  rootFile = download_assembly_summary(logger, localDir, prefix, database)

  with open(rootFile + "." + taxonomyRootId + ".idmap", "wt") as fout:
    for id in idMap:
      fout.write("%s\t%s\n" % ( id, idMap[id]))

  targetFile = extract_unique_genome(logger, rootFile, idMap, taxonomyRootId, False)
  
  localFileMap = download_assembly_genomes(logger, cacheDir, targetFile)

  categoryFile = combine_category_fasta_file(logger, maxGenomeInFile, fastaDir, prefix, targetFile, localFileMap)

  logger.info("Done.")

def prepare_database_as_whole(logger, taxonomyRootId, outputFolder, prefix, database="genbank"):
  taxonomyFile = os.path.join(outputFolder, prefix + "taxonomy.txt")
  
  prepare_taxonomy(logger, taxonomyFile)

  logger.info("Reading taxonomy from %s ..." % taxonomyFile)
  tree = TaxonomyTree()
  tree.ReadFromFile(taxonomyFile)

  idMap = tree.BuildChildIdParentNameMap(taxonomyRootId)
  name = tree.GetItem(taxonomyRootId).Name.replace(" ", "_") + "." + taxonomyRootId

  localDir = os.path.join(outputFolder, name.lower())
  if not os.path.exists(localDir):
    os.mkdir(localDir)

  cacheDir = os.path.join(localDir, "cache")
  if not os.path.exists(cacheDir):
    os.mkdir(cacheDir)

  fastaDir = os.path.join(localDir, "whole_fasta")
  if not os.path.exists(fastaDir):
    os.mkdir(fastaDir)

  rootFile = download_assembly_summary(logger, localDir, prefix, database)

  with open(rootFile + "." + taxonomyRootId + ".idmap", "wt") as fout:
    for id in idMap:
      fout.write("%s\t%s\n" % ( id, idMap[id]))

  targetFile = extract_unique_genome(logger, rootFile, idMap, taxonomyRootId, False)
  
  localFileMap = download_assembly_genomes(logger, cacheDir, targetFile)

  categoryFile = combine_fasta_file(logger, fastaDir, prefix, targetFile, localFileMap)

  logger.info("Done.")

def prepare_index(logger, categoryFile, thread, force=False, slurm=False, slurmEmail=None, slurmTemplate=None):
  bowtieIndecies = readBowtieIndexList(categoryFile)

  if not slurm:
    for bowtieIndex in bowtieIndecies:
      logger.info("Building index for %s ..." % bowtieIndex.Fasta)
      indexDone = bowtieIndex.Index + ".index.done"
      if not os.path.exists(indexDone) or force:
        subprocess.call(['bowtie-build', '-q', '-r', '--threads', str(thread), bowtieIndex.Fasta, bowtieIndex.Index])
        open(indexDone, 'wt').close()
    return

  if slurmTemplate == None:
    slurmTemplate = os.path.join(os.path.dirname(os.path.realpath(__file__)), "slurm.template")

  if not os.path.exists(slurmTemplate):
    raise Exception("Slurm template not exists: %s" % slurmTemplate)

  with open(slurmTemplate, "rt") as fin:
    slurmLines = [line.rstrip() for line in fin]

  slurmFolder = os.path.join(os.path.dirname(os.path.abspath(categoryFile)), "slurm")
  if not os.path.exists(slurmFolder):
    os.mkdir(slurmFolder)

  submitFile = os.path.join(slurmFolder, "submit.sh")
  with open(submitFile, "wt") as fout:
    for bowtieIndex in bowtieIndecies:
      indexDone = bowtieIndex.Index + ".index.done"

      if bowtieIndex.Fasta.endswith(".gz"):
        unzip_file = bowtieIndex.Fasta[:-3]
        command = "gzip -c " + bowtieIndex.Fasta + ">" + unzip_file + "\n"
      else:
        command = ""
        unzip_file = bowtieIndex.Fasta
      command = command + "bowtie-build -q -r --threads %s %s %s\necho '' > %s\n" % (thread, unzip_file, bowtieIndex.Index, indexDone )
      if bowtieIndex.Fasta.endswith(".gz"):
        command = command + "\nrm " + unzip_file + "\n"

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

def fastq_to_database(logger, fastq_file, sample_name, output_file, reads_per_file, gzipped=False):
  output_folder = os.path.dirname(os.path.abspath(output_file))
  fastaDir = os.path.join(output_folder, "fasta")
  if not os.path.exists(fastaDir):
    os.mkdir(fastaDir)

  file_index = 0
  with open(output_file, "wt") as flist:
    flist.write("BowtieIndex\tCategory\tFasta\n")
    reads_count = 0
    fout = None
    fastq_files = fastq_file.split(",")
    for fq in fastq_files:
      fin = gzip.open(fq, "rt") if fq.endswith(".gz") else open(fq, "rt")
      with fin:
        while True:
          query = fin.readline()
          if not query:
            break

          if (reads_count % reads_per_file) == 0:
            if fout != None:
              fout.close()
            file_index += 1
            bowtie_index = os.path.join(fastaDir, "%s.%d" % ( sample_name, file_index ))
            if gzipped:
              cur_file = bowtie_index + ".fasta.gz"
              fout = gzip.open(cur_file, "wt")
            else:
              cur_file = bowtie_index + ".fasta"
              fout = open(cur_file, "wt")

            flist.write("%s\t%s\t%s\n" % (bowtie_index, sample_name, cur_file ))
            logger.info("Writing to %s ..." % os.path.basename(cur_file))
          
          reads_count += 1

          sequence = fin.readline().rstrip()
          fin.readline()
          fin.readline()
          fout.write(">%d\n%s\n" %(reads_count, sequence))
    fout.close()

if __name__ == "__main__":
  logger = logging.getLogger('database')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  prepare_database(logger, '11118', "/scratch/cqs_share/references/genbank", 500, "20200331_")

