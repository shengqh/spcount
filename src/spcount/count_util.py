import argparse
import sys
import logging
import os
import gzip
import math
import subprocess
from _ctypes import ArgumentError

from .CategoryEntry import CategoryEntry
from .common_util import readFileMap
from .BowtieCountItem import BowtieCountItem, readBowtieTextFile, getQueryMap, getCategoryMap, getCategoryUniqueCountMap, assignCount, getCategoryCountMap

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
