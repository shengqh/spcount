import os
import requests
import shutil
import tarfile
from io import TextIOWrapper

from .Taxonomy import TaxonomyItem, save_taxonomy

def find_parent(item, itemMap, rank):
  result = item
  while(True):
    if result.Rank == rank:
      return result
    if result.Id == result.ParentId:
      return None
    result = itemMap[result.ParentId]

def depth(taxo, taxonomyMap):
  if taxo.Id == taxo.ParentId:
    return 0
  else:
    return depth(taxonomyMap[taxo.ParentId], taxonomyMap) + 1

def prepare_taxonomy(logger, output_file):
  #if os.path.exists(output_file):
  #  logger.info(f"File exists: {output_file}")
  #  return

  valid_ranks = ['root', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

  logger.info("Preparing taxonomy file %s ..." % output_file)

  #download taxonomy from ncbi
  url="https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
  target_file = output_file + "." + os.path.basename(url)
  if not os.path.exists(target_file):
    logger.info("Downloading %s ..." % url)
    r = requests.get(url, verify=False, stream=True)
    r.raw.decode_content = True
    with open(target_file, 'wb') as f:
      shutil.copyfileobj(r.raw, f)

  taxonomies = []

  #prepare taxonomy file with three columns: id, parentid, name
  with tarfile.open(target_file, "r:gz") as tar:
    logger.info("Reading names.dmp ...")
    with TextIOWrapper(tar.extractfile(tar.getmember('names.dmp'))) as fin:
      for line in fin:
        if "scientific name" in line:
          parts = line.split('|')
          id = int(parts[0].strip())
          name = parts[1].strip()
          taxonomies.append(TaxonomyItem(name, id, 0, ''))

    taxonomyMap = {item.Id:item for item in taxonomies}

    logger.info("Reading nodes.dmp ...")
    with TextIOWrapper(tar.extractfile(tar.getmember('nodes.dmp'))) as fin:
      for line in fin:
        parts = line.split('|')
        childId = int(parts[0].strip())
        parentId = int(parts[1].strip())
        rank = parts[2].strip()
        taxonomyMap[childId].ParentId = parentId
        taxonomyMap[childId].Rank = rank
    
    for taxo in taxonomies:
      taxo.Depth = depth(taxo, taxonomyMap)  

  logger.info("Writing result ...")
  output_ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

  for item in taxonomies:
    for rank in output_ranks:
      parentItem = find_parent(item, taxonomyMap, rank)
      item.RankMap[rank] = parentItem.Id if parentItem != None else None

  save_taxonomy(taxonomies, output_file, output_ranks)
  
  logger.info("Writing krona taxonomy db ...")
  taxonomies.sort(key=lambda x:x.Id)
  with open(os.path.join(os.path.dirname(output_file), "taxonomy.tab"), "wt") as fout:
    for taxo in taxonomies:
      if taxo.Id == taxo.ParentId:
        taxo.ParentId = 0
      fout.write(f"{taxo.Id}\t{taxo.Depth}\t{taxo.ParentId}\t{taxo.Rank}\t{taxo.Name}\n")

  #os.remove(target_file)
  logger.info(f"Result saved to {output_file} .")
