import sys

class BowtieCountItem(object):
  def __init__(self, queryName, strand, chromosome, zeroBasedOffset, category, sample="", count=1):
    self.QueryName = queryName
    self.Strand = strand
    self.Chromosome = chromosome
    self.ZeroBasedOffset = zeroBasedOffset
    self.Category = category
    self.Sample = sample
    self.Count = count
    self.Sequence = ''

def getQueryMap(bowtieItems):
  result = {}
  for item in bowtieItems:
    result.setdefault(item.QueryName, []).append(item)
  return(result)

def getCategoryMap(bowtieItems):
  result = {}
  for item in bowtieItems:
    result.setdefault(item.Category, []).append(item)
  return(result)

def getCategoryUniqueCountMap(bowtieItems):
  result = {}
  for item in bowtieItems:
    if item.Category not in result:
      result[item.Category] = 1
    else:
      result[item.Category] = result[item.Category] + 1
  return(result)

def assignCount(queryMap, bowtieItems):
  catCountMap = getCategoryUniqueCountMap(bowtieItems)
  for items in queryMap.values():
    catCounts = [catCountMap[item.Category] for item in items]
    totalCount = sum(catCounts)
    proportion = [cc / totalCount for cc in catCounts]
    for idx in range(0, len(proportion)):
      items[idx].Count = items[idx].Count * proportion[idx]

def getCategoryCountMap(bowtieItems):
  result = {}
  for item in bowtieItems:
    #print("%s\t%d" %(item.Category, item.Count))
    if item.Category not in result:
      result[item.Category] = 0
    result[item.Category] = result[item.Category] + item.Count
  return(result)

def readBowtieTextFile(fileName, sampleName):
  result = []
  with open(fileName, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      queryName = parts[0]
      result.append(BowtieCountItem(queryName, parts[1], parts[2], parts[3], parts[4], sampleName))
  return(result)
