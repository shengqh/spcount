import logging
import os
import ftplib
import gzip
import subprocess
import pandas as pd
from io import BytesIO, TextIOWrapper

# from BowtieIndex import BowtieIndexItem, readBowtieIndexList, writeBowtieIndexList
# from GenomeItem import GenomeItem, writeGenomeItems, readGenomeItems
# from Taxonomy import TaxonomyItem, TaxonomyTree

from .BowtieIndex import BowtieIndexItem, readBowtieIndexList, writeBowtieIndexList
from .GenomeItem import GenomeItem, writeGenomeItems, readGenomeItems
from .Taxonomy import TaxonomyItem, TaxonomyTree

class CategoryItem(object):
  def __init__(self, fileHandle, index, numberOfGenome):
    self.FileHandle = fileHandle
    self.Index = index
    self.NumberOfGenome = numberOfGenome

def check_file_exists(ftp, filename):
    dirname= os.path.dirname(filename)
    file_list = ftp.nlst(dirname)
    return filename in file_list

def getCategory(source):
  return(source.replace(" group", "").replace("/", "_").replace(" Bacteria", "").replace(" ", "_"))

def combine_category_fasta_file(logger, maxGenomeInFile, localDir, prefix, targetFile, localFileMap):
  categoryFile = targetFile + ".list"
  categoryFileDone = get_done_file(categoryFile)

  if os.path.exists(categoryFileDone):
    return(categoryFile)

  totalCount = len(localFileMap)

  bowtieIndecies = []

  giItems = readGenomeItems(targetFile)

  giGroup = {}
  for gi in giItems:
    category = gi.Category
    giGroup.setdefault(category, []).append(gi)

  all_categories = sorted(giGroup.keys())

  total_index = 0
  cat_index = 0
  for category in all_categories:
    cat_gis = giGroup[category]
    cat_index += 1

    subName = "%s%s_%03d" % (prefix, category, 1)
    bowtieIndex = os.path.join(localDir, subName)
    categoryFastaFile = bowtieIndex + ".fasta"
    logger.info("Init " + categoryFastaFile)

    cat = CategoryItem(open(categoryFastaFile, "wb"), 1, 1)
    bowtieIndecies.append(BowtieIndexItem(bowtieIndex, category, categoryFastaFile))

    currentCount = 0
    for gi in cat_gis:
      currentCount += 1
      total_index += 1
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
      
      logger.info(f"Merge {total_index}/{totalCount} : {cat_index}/{len(all_categories)} {category} : {currentCount}/{len(cat_gis)} : {localFnaFile} ...")
      with gzip.open(localFnaFile, 'rb') as f:
        file_content = f.read()
        categoryOut.write(file_content)
      
    cat.FileHandle.close()

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

class AssemblyGenome:
  def __init__(self, accession, taxid, organism_name, ftp_path, local_dir):
    self.accession = accession
    self.taxid = taxid
    self.organism_name = organism_name
    self.ftp_path = ftp_path
    file_prefix=os.path.basename(ftp_path)
    file_name=file_prefix + "_genomic.fna.gz"
    self.remote_fna_path=os.path.join(ftp_path, file_name)
    self.url_file = self.remote_fna_path.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
    self.url_file = self.url_file.replace("https://ftp.ncbi.nlm.nih.gov", "")
    remote_dir=os.path.dirname(self.url_file)
    sub_folder = os.path.join(local_dir, remote_dir[1:])
    if not os.path.exists(sub_folder):
      os.makedirs(sub_folder)
    self.local_fna_path=os.path.join(sub_folder, file_name)
    self.local_done_path=self.local_fna_path + ".done"

def download_assembly_summary(logger, output_file, database="refseq"):
  if os.path.exists(output_file):
    logger.info(f"File exists: {output_file}")
    return

  summaryFile =  "/genomes/%s/assembly_summary_%s.txt" % (database, database)

  logger.info("Downloading %s ..." % summaryFile)
  with open_ftp() as ftp:
    with open(output_file, "wb") as f:
      ftp.retrbinary("RETR " + summaryFile, f.write)

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
  logger.info("Reading %s ... " % targetFile)

  giItems = readGenomeItems(targetFile)
  totalCount = len(giItems)

  logger.info("Checking %d genomes in cache folder %s ..." % (totalCount, genomeRootDir))
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

        remoteFile=fnaFile.replace("https://ftp.ncbi.nlm.nih.gov", "")

        logger.info("Downloading %d/%d: %s to %s ..." % (currentCount, totalCount, remoteFile, localFile))
        for retry in [1,2,3]:
          #time.sleep(0.1)
          try:
            with open(localFile, "wb") as f:
              ftp.retrbinary("RETR " + remoteFile, f.write, 1024)
              open(localDoneFile, 'wt').close()
          except:
            logger.error(f"ERROR: downloading {remoteFile} retry {retry} failed.")
            ftp.close()
            ftp = open_ftp()

  return(localFileMap)

class AssemblyGenome:
  def __init__(self, accession, taxid, organism_name, ftp_path, local_dir):
    self.accession = accession
    self.taxid = taxid
    self.organism_name = organism_name
    self.ftp_path = ftp_path
    file_prefix=os.path.basename(ftp_path)
    file_name=file_prefix + "_genomic.fna.gz"
    self.remote_fna_path=os.path.join(ftp_path, file_name)
    self.url_file = self.remote_fna_path.replace("ftp://ftp.ncbi.nlm.nih.gov", "")
    self.url_file = self.url_file.replace("https://ftp.ncbi.nlm.nih.gov", "")
    remote_dir=os.path.dirname(self.url_file)
    sub_folder = os.path.join(local_dir, remote_dir[1:])
    if not os.path.exists(sub_folder):
      os.makedirs(sub_folder)
    self.local_fna_path=os.path.join(sub_folder, file_name)
    self.local_done_path=self.local_fna_path + ".done"

def prepare_segment_database(logger, taxonomyFile, assemblySummaryFile, taxonomyRootId, outputFolder, prefix, genomeNumberPerFile=500, referenceAndRepresentativeOnly=True):
  # taxonomyFile = '/data/cqs/references/bacteria/20220406_taxonomy.txt'
  # assemblySummaryFile = '/data/cqs/references/bacteria/20220406_assembly_summary_refseq.txt'
  # taxonomyRootId = 2 
  # outputFolder = '/data1/shengq2/references/spcount'
  # prefix = '20220406'

  ref_categories = ['reference genome', 'representative genome']

  logger.info(f"Reading taxonomy from {taxonomyFile} ...")
  taxonomy=pd.read_csv(taxonomyFile, sep="\t", index_col=0)

  if not taxonomyRootId in taxonomy.index:
    raise Exception(f'Cannot find taxonomy id {taxonomyRootId} in {taxonomyFile}')

  root=taxonomy.loc[taxonomyRootId]
  root_taxonomy=taxonomy[taxonomy[root.Rank]==taxonomyRootId]
  
  localDir = outputFolder

  cacheDir = os.path.join(localDir, "cache")
  if not os.path.exists(cacheDir):
    os.mkdir(cacheDir)

  fastaDir = os.path.join(localDir, "fasta")
  if not os.path.exists(fastaDir):
    os.mkdir(fastaDir)

  logger.info(f"Reading assembly summary from {assemblySummaryFile} ...")
  assembly=pd.read_csv(assemblySummaryFile, sep="\t", header=1, index_col=0)

  root_assembly=assembly[assembly.taxid.isin(root_taxonomy.index)]
  if referenceAndRepresentativeOnly:
    root_assembly=root_assembly.loc[root_assembly.refseq_category.isin(ref_categories)]

  genomes = []
  for row in root_assembly.itertuples():
    genome = AssemblyGenome(row.Index, row.taxid, row.organism_name, row.ftp_path, cacheDir )
    genomes.append(genome)
  logger.info(f"Total {len(genomes)} genomes ...")

  for rep in [1,2,3]:
    #cache all genome
    missing_genomes =[g for g in genomes if not os.path.exists(g.local_done_path)] 
    missing_count = len(missing_genomes)
    if missing_count > 0:
      logger.info(f"Downloading {missing_count} out of {len(genomes)} genomes ...")
      with open_ftp() as ftp:
        currentCount = 0
        for genome in missing_genomes:
          currentCount = currentCount + 1

          localFile = genome.local_fna_path
          remoteFile = genome.url_file

          logger.info("Downloading %d/%d: %s ..." % (currentCount, missing_count, os.path.basename(localFile)))
          for retry in [1,2,3]:
            #time.sleep(0.1)
            try:
              with open(localFile, "wb") as f:
                ftp.retrbinary("RETR " + remoteFile, f.write, 1024)
                open(genome.local_done_path, 'wt').close()
                break
            except:
              logger.error(f"Error: downloading {remoteFile} retry {retry} failed.")
              ftp.close()
              ftp = open_ftp()

  #download gtf file
  logger.error("Downloading gtf files.")
  with open_ftp() as ftp:
    currentCount = 0
    for genome in genomes:
      currentCount = currentCount + 1
      localFile = genome.local_fna_path
      remoteFile = genome.url_file

      if os.path.exists(genome.local_done_path):
        remoteGtfFile = remoteFile.replace("_genomic.fna.gz", "_genomic.gtf.gz")
        if check_file_exists(ftp, remoteGtfFile):
          localGtfFile = localFile.replace("_genomic.fna.gz", "_genomic.gtf.gz")
          if not os.path.exists(localGtfFile):
            logger.info(f"Downloading {currentCount}/{len(genomes)}: {os.path.basename(localGtfFile)} ...")
            with open(localGtfFile, "wb") as f:
              ftp.retrbinary("RETR " + remoteGtfFile, f.write, 1024)
        else:
          logger.error(f"Remote file {remoteGtfFile} not exists.")
          #raise Exception(f"Remote file {remoteGtfFile} not exists.")

  missing_genomes =[g for g in genomes if not os.path.exists(g.local_done_path)] 
  missing_count = len(missing_genomes)

  if len(missing_genomes) > 0:
    logger.error("After multiple tries, there are still some genomes failed.")
    for miss in missing_genomes:
      logger.error(f"  {miss.remote_file_path}")
  else:
    logger.info("Output fasta file ...")
    output_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    bowtieIndecies = []
    with open(os.path.join(outputFolder, prefix + ".taxonomy.txt"), "wt") as ftx:
      ftx.write("chrom\taccession\tscientific_name\ttaxid\trank\t%s\n" % "\t".join(output_ranks))
      findex = 0
      fFasta = None
      gindex = 0
      for genome in genomes:
        gtax = taxonomy.loc[genome.taxid]
        if gindex % genomeNumberPerFile == 0:
          if fFasta != None:
            fFasta.close()
          findex += 1
          bowtie_index = os.path.join(fastaDir, "%s.%03d" % (prefix, findex))
          cur_file = bowtie_index + ".fa"
          bowtieIndecies.append(BowtieIndexItem(bowtie_index, root.ScientificName, cur_file))
          logger.info(f"Writing to {gindex+1}/{len(genomes)}: {cur_file} ...")
          fFasta = open(cur_file, "wb")
        gindex += 1
        with gzip.open(genome.local_fna_path, 'rb') as f:
          file_content = f.read()
          fFasta.write(file_content)
          with TextIOWrapper(BytesIO(file_content)) as fin:
            for line in fin:
              if line.startswith(">"):
                parts = line.rstrip().split(' ')
                chrom = parts[0][1:]
                ftx.write("%s\t%s\t%s\t%s\t%s" % (chrom, genome.accession, genome.organism_name, genome.taxid, gtax['Rank']))
                for rank in output_ranks:
                  taxonomyId = gtax[rank]
                  if pd.isnull(taxonomyId):
                    ftx.write("\tUnclassified")
                  else:
                    rtex=taxonomy.loc[taxonomyId]
                    ftx.write(f"\t{rtex.ScientificName}")
                ftx.write("\n")
        fFasta.close() 
    writeBowtieIndexList(os.path.join(outputFolder, prefix + ".index.txt"), bowtieIndecies)

    gtf_file = os.path.join(fastaDir, prefix + ".gtf")
    with open(gtf_file, "wb") as fgtf:
      with open(gtf_file + ".missing", "wt") as fgtfmiss:
        for genome in genomes:
          local_gtf_file = genome.local_fna_path.replace("_genomic.fna.gz", "_genomic.gtf.gz")
          if(os.path.exists(local_gtf_file)):
            with gzip.open(local_gtf_file, 'rb') as f:
              file_content = f.read()
              fgtf.write(file_content)
          else:
            fgtfmiss.write(f"{genome.local_fna_path}\n")

  logger.info("Done.")

def prepare_index(logger, categoryFile, thread, force=False, slurmTemplate=None):
  bowtieIndecies = readBowtieIndexList(categoryFile)

  if slurmTemplate == None:
    for bowtieIndex in bowtieIndecies:
      logger.info("Building index for %s ..." % bowtieIndex.Fasta)
      indexDone = bowtieIndex.Index + ".index.done"
      if not os.path.exists(indexDone) or force:
        subprocess.call(['bowtie-build', '-q', '-r', '--threads', str(thread), bowtieIndex.Fasta, bowtieIndex.Index])
        open(indexDone, 'wt').close()
    return

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

      command = command + """
rm -f %s.index.failed

bowtie-build -q -r --threads %s %s %s

status=$?
if [[ $status -ne 0 ]]; then
  touch %s.index.failed
else
  touch %s
fi 
""" % (bowtieIndex.Index, thread, unzip_file, bowtieIndex.Index, bowtieIndex.Index, indexDone )

      if bowtieIndex.Fasta.endswith(".gz"):
        command = command + "\nrm " + unzip_file + "\n"

      slurmFile = os.path.join(slurmFolder, os.path.basename(bowtieIndex.Index + ".slurm"))
      with open(slurmFile, "wt") as fslurm:
        for line in slurmLines:
          if "__THREAD__" in line:
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
  #prepare_database(logger, '11118', "/data/cqs/references/spcount", 500, "20211111_")
  prepare_database(logger, '2', "/data/cqs/references/spcount", 500, "20211111_")

