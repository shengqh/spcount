import gzip
import re
import pickle

from .CategoryEntry import CategoryEntry
from .common_util import readFileMap
from .BowtieCountItem import BowtieCountItem, readBowtieTextFile, getQueryMap, assignCount
from .Species import Species
from .Query import Query

def removeSubset(logger, catMap):
  catItems = [v for v in catMap.values()]
  catItems.sort(key=lambda c:len(c.Items), reverse=True)
  for ci in catItems:
    ci.QueryNames = set(bi.QueryName for bi in ci.Items)

  logger.info("Removing subset ...")
  bRemoved = False
  for si in range(len(catItems) - 1, 1, -1):
    sQueryNames = catItems[si].QueryNames
    for hi in range(si-1, 0, -1):
      hQueryNames = catItems[hi].QueryNames
      if sQueryNames.issubset(hQueryNames):
        logger.info("  Remove %s as subset of %s" % (catItems[si].Name, catItems[hi].Name))
        bRemoved = True
        del catItems[si]
        break

  if not bRemoved:
    return(catMap)

  result = {cat.Name:cat for cat in catItems }
  return(result)

def fillQueryCount(logger, queryMap, countFile):
  logger.info("Filling count from %s ..." % countFile)
  with open(countFile, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      queryName = parts[0]
      
      if not queryName in queryMap:
        continue
      
      count = float(parts[1])
      sequence = parts[2]

      items = queryMap[queryName]
      for item in items:
        item.Count = count
        item.Sequence = sequence  

def bowtie_count(logger, input_list_file, output_file, count_file, species_file, species_column='species'):
  logger.info(f"reading count file {count_file}")
  count_map={}
  with open(count_file, "rt") as fin:
    fin.readline()
    for line in fin:
      parts=line.split('\t')
      count_map[parts[0]] = parts[1]

  logger.info(f"reading species file {species_file}")
  species_map={}
  with open(species_file, "rt") as fin:
    line = fin.readline()
    headers = line.rstrip().split('\t')
    species_index = headers.index(species_column)
    for line in fin:
      parts = line.rstrip().split('\t')
      species_map[parts[0]] = parts[species_index]

  read_map = {}
  with open(input_list_file, "rt") as fl:
    for line in fl:
      parts = re.split('\s+', line.rstrip())
      bowtie_file = parts[0]
      logger.info(f"parsing {bowtie_file}")
      with gzip.open(bowtie_file, "rt") as fin:
        for bl in fin:
          bparts = bl.split('\t')
          query = bparts[0].split(' ')[0]
          species = species_map[bparts[2]]
          read_map.setdefault(query,set()).add(species)

  logger.info(f"merge all bowtie result ...")
  all_queries = [[query, int(count_map[query]), ",".join(read_map[query])] for query in read_map.keys()]   
  all_queries.sort(key=lambda x:x[1], reverse=True)

  logger.info(f"output to {output_file} ...")
  with gzip.open(output_file, "wt") as fout:
    fout.write("read\tcount\tspecies\n")
    for query in all_queries:
      fout.write(f"{query[0]}\t{query[1]}\t{query[2]}\n")

  logger.info("done")

def count_table(logger, input_list_file, output_prefix, species_file, species_column='species'):

  logger.info("reading species taxonomy map from " + species_file + "...")
  species_taxonomy_map = {}
  with open(species_file, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      if parts[12] not in species_taxonomy_map:
        species_taxonomy_map[parts[12]] = {
          'genus':parts[11],
          'family':parts[10],
          'order':parts[9],
          'class':parts[8],
          'phylum':parts[7],
          'superkingdom':parts[5],
      }

  query_list = []
  species_map = {}
  samples=[]
  with open(input_list_file, "rt") as fl:
    for line in fl:
      parts = line.rstrip().split('\t')
      count_file = parts[0]
      sample=parts[1]
      samples.append(sample)
      logger.info(f"parsing {count_file}")
      with gzip.open(count_file, "rt") as fin:
        fin.readline()
        for bl in fin:
          bparts = bl.split('\t')
          query_name = bparts[0].split(' ')[0]
          count = int(bparts[1])
          species_list = bparts[2].rstrip().split(',')
          query = Query(query_name, count, species_list)
          query_list.append(query)
          for species in species_list:
            if species not in species_map:
              species_map[species] = Species(species)
            species_map[species].add_auery(sample, query)

  species_list = [v for v in species_map.values()]
  species_map.clear()

  for sv in species_list:
    sv.sum_query_count()

  species_list.sort(key=lambda x:x.query_count, reverse=True)

  # myobj={ "species_taxonomy_map":species_taxonomy_map, 
  #         "species_list": species_list,
  #         "samples": samples,
  #         "query_list": query_list }

  # with open(output_prefix + '.pickle', 'wb') as handle:
  #   pickle.dump(myobj, handle, protocol=pickle.HIGHEST_PROTOCOL)

  logger.info("merge identical")
  for sv in species_list:
    sv.is_identical = False
    sv.identical_species = []

  for i1 in range(0, len(species_list)-1):
    if species_list[i1].is_identical:
      continue
    if i1 % 100 == 0:
      logger.info(f"checking identical: {i1+1} / {len(species_list)}")
    for i2 in range(i1+1, len(species_list)):
      if species_list[i2].is_identical:
        continue
      if species_list[i1].query_count != species_list[i2].query_count:
        break
      if species_list[i1].queries_set == species_list[i2].queries_set:
        species_list[i2].is_identical = True
        species_list[i1].identical_species.append(species_list[i2])

  old_len = len(species_list)
  new_len = len([sv for sv in species_list if not sv.is_identical])
  logger.info(f"{old_len - new_len} identical species were found")

  for i1 in range(0, len(species_list)-1):
    if species_list[i1].is_subset or species_list[i1].is_identical:
      continue
    if i1 % 100 == 0:
      logger.info(f"checking subset: {i1+1} / {len(species_list)}")
    for i2 in range(i1+1, len(species_list)):
      if species_list[i2].is_subset or species_list[i1].is_identical:
        continue
      if species_list[i1].contains(species_list[i2]):
        species_list[i2].is_subset = True

  #remove subset species from query 
  for species in species_list:
    if species.is_subset:
      all_names = [species.name] + species.identical_species
      for qlist in species.queries.values():
        for q in qlist:
          q.remove_species(all_names)

  old_len = len(species_list)
  species_list = [sv for sv in species_list if not sv.is_subset]
  new_len = len(species_list)
  logger.info(f"{old_len - new_len} subset were removed")

  with open(output_prefix + ".species.query.count", "wt") as fout:
    fout.write("Feature\t" + "\t".join(samples) + "\n")
    for species in species_list:
      if species.is_identical:
        continue
      species_name = species.name
      if len(species.identical_species) > 0:
        species_name = species_name + "," + ",".join([s.name for s in species.identical_species])
      countstr = "\t".join(str(species.sample_query_count[sample]) if sample in species.sample_query_count else "0" for sample in samples)
      fout.write(f"{species_name}\t{countstr}\n")

  for query in query_list:
    query.estimate_count()

  for species in species_list:
    species.sum_estimated_count()

  with open(output_prefix + ".species.estimated.count", "wt") as fout:
    fout.write("Feature\t" + "\t".join(samples) + "\n")
    for species in species_list:
      if species.is_identical:
        continue

      species_name = species.name
      if len(species.identical_species) > 0:
        species_name = species_name + "," + ",".join([s.name for s in species.identical_species])
      countstr = "\t".join("{:.2f}".format(species.sample_estimated_count[sample]) if sample in species.sample_estimated_count else "0" for sample in samples)
      fout.write(f"{species_name}\t{countstr}\n")

  levels = [ 'genus', 'family', 'order', 'class', 'phylum']
  for level in levels:
    logger.info("output " + level)
    cat_map = {}
    for species in species_list:
      cat_name = species_taxonomy_map[species.name][level]
      #print(species.name + ": " + cat_name)
      if cat_name not in cat_map:
        cat_map[cat_name] = Species(cat_name)
      cat_map[cat_name].identical_species.append(species)
    
    cats = [cat for cat in cat_map.values()]
    for cat in cats:
      cat.sample_query_count = {}
      for sample in samples:
        squeries_set = set()
        for species in cat.identical_species:
          if sample in species.queries:
            for query in species.queries[sample]:
              squeries_set.add(query)
        scount = sum([query.count for query in squeries_set])
        cat.sample_query_count[sample] = scount
      cat.query_count = sum([v for v in cat.sample_query_count.values()])

      cat.estimated_count = sum([species.estimated_count for species in cat.identical_species])
      cat.sample_estimated_count = {sample: sum([species.sample_estimated_count[sample] if sample in species.sample_estimated_count else 0 for species in cat.identical_species]) for sample in samples}
      
    cats.sort(key=lambda x:x.query_count, reverse=True)

    with open(output_prefix + "." + level + ".query.count", "wt") as fout:
      fout.write("Feature\t" + "\t".join(samples) + "\n")
      for species in cats:
        species_name = species.name
        countstr = "\t".join(str(species.sample_query_count[sample]) if sample in species.sample_query_count else "0" for sample in samples)
        fout.write(f"{species_name}\t{countstr}\n")

    with open(output_prefix + "." + level + ".estimated.count", "wt") as fout:
      fout.write("Feature\t" + "\t".join(samples) + "\n")
      for species in cats:
        species_name = species.name
        countstr = "\t".join("{:.2f}".format(species.sample_estimated_count[sample]) if sample in species.sample_estimated_count else "0" for sample in samples)
        fout.write(f"{species_name}\t{countstr}\n")

  logger.info("done")

def count(logger, inputListFile, outputFile, countListFile, category_name=None):
  logger.info("Start count ...")

  bowtieFileMap = readFileMap(inputListFile)

  #print(bowtieFileMap)

  countFileMap = readFileMap(countListFile) if countListFile != None else {}

  # icount = 0
  finalItems = []
  for sample in bowtieFileMap.keys():
    bowtieFile = bowtieFileMap[sample]

    #for debug
    #if not os.path.exists(bowtieFile):
    #  continue

    logger.info("Processing %s ..." % bowtieFile)

    bowtieItems = readBowtieTextFile(bowtieFile, sample)
     
    queryMap = getQueryMap(bowtieItems)

    if sample in countFileMap:
      countFile = countFileMap[sample]
      fillQueryCount(logger, queryMap, countFile)

    assignCount(queryMap, bowtieItems)

    finalItems.extend(bowtieItems)

  samples = sorted(list(set([bi.Sample for bi in finalItems])))

  if category_name != None:
    for bi in finalItems:
      bi.Category = category_name

  categorySet = list(set([bi.Category for bi in finalItems]))

  finalMap = {sample:{category:0 for category in categorySet} for sample in samples}
  for bi in finalItems:
    finalMap[bi.Sample][bi.Category] += bi.Count
  
  catCount = [[cat, sum(v[cat] for v in finalMap.values())] for cat in categorySet]
  catCount.sort(key=lambda r:-r[1])

  categories = sorted(list(cat[0] for cat in catCount))
  
  logger.info("Writing to %s" % outputFile)
  with open(outputFile, "wt") as fout:
    fout.write("Category\t%s\n" % "\t".join(samples))
    for catName in categories:
      fout.write(catName)
      for sample in samples:
        catMap = finalMap[sample]
        if catName in catMap:
          fout.write("\t%.2lf" % finalMap[sample][catName])
        else:
          fout.write("\t0")
      fout.write("\n")

  logger.info("done")

def count_old(logger, inputListFile, outputFile, countListFile):
  logger.info("Start count ...")

  bowtieFileMap = readFileMap(inputListFile)
  countFileMap = readFileMap(countListFile) if countListFile != None else {}

  # icount = 0
  finalMap = {}
  for sample in bowtieFileMap.keys():
    bowtieFile = bowtieFileMap[sample]

    #for debug
    #if not os.path.exists(bowtieFile):
    #  continue

    logger.info("Processing %s ..." % bowtieFile)

    bowtieItems = []
    with open(bowtieFile, "rt") as fin:
      for line in fin:
        parts = line.rstrip().split('\t')
        queryName = parts[0]
        bowtieItems.append(BowtieCountItem(queryName, parts[1], parts[2], parts[3], parts[4]))

    catMap = {}
    for item in bowtieItems:
      catMap.setdefault(item.Category, CategoryEntry(item.Category)).Items.append(item)

    #catMap = removeSubset(logger, catMap)

    catItems = [v for v in catMap.values()]
    catItems.sort(key=lambda c:len(c.Items), reverse=True)
    for ci in catItems:
      logger.info("%s : %d" % (ci.Name, len(ci.Items)))

    queryMap = {}
    for cat in catMap.values():
      for item in cat.Items:
        queryMap.setdefault(item.QueryName, []).append(item)

    if sample in countFileMap:
      countFile = countFileMap[sample]
      fillQueryCount(logger, queryMap, countFile)

    for itemList in queryMap.values():
      count = itemList[0].Count
      subCount = count / len(itemList)
      for item in itemList:
        item.Count = subCount

    for ci in catMap.values():
      ci.TotalCount = sum([bi.Count for bi in ci.Items])

    finalMap[sample] = catMap
    # icount = icount + 1
    # if icount ==2 :
    #   break

  samples = sorted([s for s in finalMap.keys()])

  categorySet = set()
  for catMap in finalMap.values():
    categorySet = categorySet.union([cat for cat in catMap.keys()])
  
  catCount = [[cat, sum(v[cat].TotalCount for v in finalMap.values() if cat in v)] for cat in categorySet]
  catCount.sort(key=lambda r:-r[1])

  categories = list(cat[0] for cat in catCount)
  
  logger.info("Writing to %s" % outputFile)
  with open(outputFile, "wt") as fout:
    fout.write("Category\t%s\n" % "\t".join(samples))
    for catName in categories:
      fout.write(catName)
      for sample in samples:
        catMap = finalMap[sample]
        if catName in catMap:
          fout.write("\t%.2lf" % finalMap[sample][catName].TotalCount)
        else:
          fout.write("\t0")
      fout.write("\n")

  logger.info("done")
