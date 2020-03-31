import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from _ctypes import ArgumentError

def readFileMap(fileName):
  result = {}
  with open(fileName, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      result[parts[1]] = parts[0]
  return(result)

class BowtieCountItem(object):
  def __init__(self, queryName, strand, chromosome, zeroBasedOffset, category, sample, count=1):
    self.QueryName = queryName
    self.Strand = strand
    self.Chromosome = chromosome
    self.ZeroBasedOffset = zeroBasedOffset
    self.Category = category
    self.Sample = sample
    self.Count = count
    self.Sequence = ''

class CategoryEntry(object):
  def __init__(self, name):
    self.Name = name
    self.TotalCount = 0
    self.Items = []
    self.QueryNames = set()

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

def count(logger, inputListFile, outputFile, countListFile):
  logger.info("Start count ...")

  bowtieFileMap = readFileMap(inputListFile)
  countFileMap = readFileMap(countListFile) if countListFile != None else {}

  # icount = 0
  finalItems = []
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
        bowtieItems.append(BowtieCountItem(queryName, parts[1], parts[2], parts[3], parts[4], sample))

    queryMap = getQueryMap(bowtieItems)

    if sample in countFileMap:
      countFile = countFileMap[sample]
      logger.info("filling count from %s ..." % countFile)
      fillQueryCount(logger, queryMap, countFile)

    assignCount(logger, queryMap, bowtieItems)

    finalItems.extends(bowtieItems)

  samples = sorted(list(set([bi.Sample for bi in finalItems])))

  categorySet = list(set([bi.Category for bi in finalItems]))
  
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
