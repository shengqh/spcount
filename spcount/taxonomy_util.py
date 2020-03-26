import gzip
import os
import sys
import requests
import shutil
import tarfile
from io import TextIOWrapper
from .Taxonomy import TaxonomyItem

def prepare_taxonomy(logger, outputFile, force=False):
  if not os.path.exists(outputFile) or force:
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

  logger.info("done")

