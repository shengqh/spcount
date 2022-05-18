import gzip
import re
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import OrderedDict

from .CategoryEntry import CategoryEntry
from .common_util import readFileMap
from .BowtieCountItem import BowtieCountItem, readBowtieTextFile, getQueryMap, assignCount
from .Species import Species
from .Query import Query
from .Sequence import Sequence

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
      parts=line.rstrip().split('\t')
      count_map[parts[0]] = [parts[1], parts[2]]

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
  all_queries = [[query, int(count_map[query][0]), count_map[query][1], ",".join(read_map[query])] for query in read_map.keys()]   
  all_queries.sort(key=lambda x:x[1], reverse=True)

  logger.info(f"output to {output_file} ...")
  with gzip.open(output_file, "wt") as fout:
    fout.write("read\tcount\tsequence\tspecies\n")
    for query in all_queries:
      fout.write(f"{query[0]}\t{query[1]}\t{query[2]}\t{query[3]}\n")

  logger.info("done")

def read_species_taxonomy_map(species_file):
  species_taxonomy_map = {}
  with open(species_file, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      if parts[12] not in species_taxonomy_map:
        species_taxonomy_map[parts[12]] = {
          'species':parts[12],
          'genus':parts[11],
          'family':parts[10],
          'order':parts[9],
          'class':parts[8],
          'phylum':parts[7],
          'superkingdom':parts[5],
      }
  return(species_taxonomy_map)

def read_file_map(input_list_file):
  file_map = OrderedDict()
  with open(input_list_file, "rt") as fl:
    for line in fl:
      parts = line.rstrip().split('\t')
      file_map[parts[1]]=parts[0]
  return(file_map)

def read_sequence_list(logger, file_map, debug_mode=False):
  sequence_map = {}
  for sample, count_file in file_map.items():
    if logger != None:
      logger.info(f"parsing {count_file} for sequence ...")

    with gzip.open(count_file, "rt") as fin:
      fin.readline()
      bcount = 0
      for bl in fin:
        bparts = bl.split('\t')
        count = int(bparts[1])
        seq = bparts[2]
        sequence = sequence_map.get(seq)
        if sequence == None:
          sequence = Sequence(seq)
          sequence_map[seq] = sequence
        sequence.add_sample_count(sample, count)

        bcount += 1
        if debug_mode and bcount == 10000:
          break
  result = list(sequence_map.values())
  result.sort(key=lambda x:x.query_count, reverse=True)
  return(result)

def read_query_list(logger, file_map, debug_mode=False):
  query_list = []
  for sample, count_file in file_map.items():
    if logger != None:
      logger.info(f"parsing {count_file} for query")

    with gzip.open(count_file, "rt") as fin:
      fin.readline()
      bcount = 0
      for bl in fin:
        bparts = bl.split('\t')
        query_name = bparts[0].split(' ')[0]
        count = int(bparts[1])
        species_list = bparts[3].rstrip().split(',')
        query = Query(sample, query_name, count, species_list)
        query_list.append(query)
        bcount += 1
        if debug_mode and bcount == 10000:
          break
  return(query_list)

def build_species_list(query_list):
  species_map = {}
  for query in query_list:
    for species in query.species_list:
      if species not in species_map:
        species_map[species] = Species(species)
      species_map[species].add_query(query)
  result = list(species_map.values())
  for species in result:
    species.sum_query_count()
  result.sort(key=lambda x:x.query_count, reverse=True)
  return(result)

# aggregate each query to unique node. The result contains nodes from different ranks.
def build_aggregate_node_list(query_list, species_taxonomy_map, ranks, aggregate_rate):
  agg_species_list = {}
  for q in query_list:
    rank, rank_name = q.get_unique_rank(species_taxonomy_map, ranks, aggregate_rate)
    species = agg_species_list.get(rank_name)
    if species == None:
      species = Species(rank_name, rank)
      agg_species_list[rank_name] = species
    species.add_query(q)

  result = list(agg_species_list.values())
  for species in result:
    species.sum_query_count()
  result.sort(key=lambda x:x.query_count, reverse=True)
  return(result)

# aggregate each query to specific rank. The query doesn't aggregated to one node will be assigned to None.
def build_aggregate_rank_list(query_list, species_taxonomy_map, rank, aggregate_rate):
  rank_map = {}
  for q in query_list:
    rank_name = q.aggregate_to_rank(species_taxonomy_map, rank, aggregate_rate)
    rank_obj = rank_map.get(rank_name)
    if rank_obj == None:
      rank_obj = Species(rank_name, rank)
      rank_map[rank_name] = rank_obj
    rank_obj.add_query(q)
  result = list(rank_map.values())
  for rank_obj in result:
    rank_obj.sum_query_count()
  result.sort(key=lambda x:x.query_count, reverse=True)
  return(result)

def output_rank_list(output_file, rank_list, samples, with_tax_id=False):
  with open(output_file, "wt") as fout:
    if with_tax_id:
      fout.write("Rank_name\tTaxonomyId\tRank\t" + "\t".join(samples) + "\n")
    else:
      fout.write("Rank_name\tRank\t" + "\t".join(samples) + "\n")
    for rank_obj in rank_list:
      sample_count_str = rank_obj.get_query_count_str(samples)
      if with_tax_id:
        fout.write(f"{rank_obj.name}\t{rank_obj.taxid}\t{rank_obj.rank}\t{sample_count_str}\n")
      else:
        fout.write(f"{rank_obj.name}\t{rank_obj.rank}\t{sample_count_str}\n")
  
def count_table(logger, input_list_file, output_prefix, taxonomy_file, species_file, species_column='species', aggregate_rate=0.95, debug_mode=False):
  logger.info(f"Reading taxonomy from {taxonomy_file} ...")
  taxonomy=pd.read_csv(taxonomy_file, sep="\t")
  taxonomy_name_id_map=dict(zip(taxonomy.ScientificName, taxonomy.Id))
  taxonomy_name_id_map['AmbiguousRanks'] = -1
  taxonomy_name_id_map['Unclassified'] = -1

  logger.info("reading species taxonomy map from " + species_file + "...")
  species_taxonomy_map = read_species_taxonomy_map(species_file)

  file_map = read_file_map(input_list_file)

  samples=list(file_map.keys())

  #in order to save memory, we handle the sequence first. then in query mode, we don't need to store sequence anymore.
  logger.info("building sequence list ...")
  sequence_list = read_sequence_list(logger, file_map, debug_mode)
  with open(output_prefix + ".read.count", "wt") as fout:
    fout.write("Sequence\t" + "\t".join(samples) + "\n")
    for sequence in sequence_list:
      sample_count_str = sequence.get_sample_count_str(samples)
      fout.write(f"{sequence.seq}\t{sample_count_str}\n")
  sequence_list = None

  logger.info("building query list ...")
  query_list = read_query_list(logger, file_map, debug_mode)

  logger.info("building species list from query list ...")
  species_list = build_species_list(query_list)

  for i1 in range(0, len(species_list)-1):
    if species_list[i1].is_subset:
      continue
    if i1 % 100 == 0:
      logger.info(f"checking subset: {i1+1} / {len(species_list)}")
    for i2 in range(i1+1, len(species_list)):
      if species_list[i2].is_subset:
        continue
      if species_list[i1].contains(species_list[i2]):
        species_list[i2].is_subset = True

  logger.info("remove subset species from query")
  for species in species_list:
    if species.is_subset:
      for qlist in species.queries.values():
        for q in qlist:
          q.species_list.remove(species.name)

  old_len = len(species_list)
  species_list = [sv for sv in species_list if not sv.is_subset]
  new_len = len(species_list)
  logger.info(f"{old_len - new_len} subset were removed")

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
  # if debug_mode:
  #   myobj={ "species_taxonomy_map":species_taxonomy_map, 
  #           "species_list": species_list,
  #           "samples": samples,
  #           "query_list": query_list }
  #   with open(output_prefix + '.pickle', 'wb') as handle:
  #     pickle.dump(myobj, handle, protocol=pickle.HIGHEST_PROTOCOL)

  logger.info(f"output aggregated count of rank species ...")
  rank_list = build_aggregate_rank_list(query_list, species_taxonomy_map, "species", aggregate_rate)
  output_rank_list(output_prefix + f".species.aggregated.count", rank_list, samples, with_tax_id=False)
  rank_list = None

  logger.info(f"output query count of rank species ...")
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

  logger.info(f"output estimated count of rank species ...")
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
    logger.info(f"output aggregated count of rank {level} ...")
    rank_list = build_aggregate_rank_list(query_list, species_taxonomy_map, level, aggregate_rate)
    output_rank_list(output_prefix + f".{level}.aggregated.count", rank_list, samples, with_tax_id=False)
    rank_list = None

    logger.info(f"output query/estimated count of rank {level} ...")
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
        countstr = species.get_query_count_str(samples)
        fout.write(f"{species_name}\t{countstr}\n")

    with open(output_prefix + "." + level + ".estimated.count", "wt") as fout:
      fout.write("Feature\t" + "\t".join(samples) + "\n")
      for species in cats:
        species_name = species.name
        countstr = species.get_estimated_count_str(samples)
        fout.write(f"{species_name}\t{countstr}\n")

  logger.info("output aggregated node ...")
  ranks=[ 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
  rank_list = build_aggregate_node_list(query_list, species_taxonomy_map, ranks, aggregate_rate)
  for rank in rank_list:
    rank.taxid = taxonomy_name_id_map[rank.name]
  output_rank_list(output_prefix + ".tree.count", rank_list, samples, with_tax_id=True)

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
