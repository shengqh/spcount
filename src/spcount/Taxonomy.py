import unittest
import requests
import os
import shutil
import tarfile
from collections import deque

class TaxonomyItem(object):
  def __init__(self, name, id, parentid, rank):
    self.Name = name
    self.Id = id
    self.ParentId = parentid
    self.Rank = rank
    self.RankMap = {}

  def __str__(self):
    return (f"{self.Name},{self.Id},{self.Rank}")

def save_taxonomy(taxonomies, output_file, output_ranks):
  with open(output_file, "wt") as fout:
    fout.write("Id\tParentId\tScientificName\tRank\t%s\n" % ("\t".join(output_ranks)))
    for item in taxonomies:
      fout.write(f"{item.Id}\t{item.ParentId}\t{item.Name}\t{item.Rank}")
      for rank in output_ranks:
        rank_value = item.RankMap[rank]
        if rank_value == None:
          fout.write("\t")
        else:
          fout.write(f"\t{rank_value}")
      fout.write("\n")

def load_taxonomy_map(input_file, numberOfEntry=-1):
  items = []
  with open(input_file, "rt") as fin:
    header = fin.readline()
    header_parts = header.rstrip().split('\t')

    count = 0
    for line in fin:
      parts = line.split('\t')
      parts[len(parts)-1] = parts[len(parts)-1].rstrip()
      item = TaxonomyItem(parts[2], int(parts[0]), int(parts[1]), parts[3])
      for idx in range(4, len(header_parts)):
        rank_name = header_parts[idx]
        rank_value = parts[idx]
        item.RankMap[rank_name] = int(rank_value) if rank_value != '' else None
      items.append(item)

      count += 1
      if numberOfEntry > 0 and count > numberOfEntry:
        break
    
  idItemMap = {item.Id:item for item in items}

  return idItemMap

class TaxonomyTree(object):
  def __init__(self):
    self.IdItemMap = {}
    self.ChildToParentMap = {}
    self.ParentChildTree = {}

  def ReadFromFile(self, taxonomyFile, numberOfEntry=-1):
    self.ChildToParentMap.clear()
    self.ParentChildTree.clear()

    self.IdItemMap = load_taxonomy_map(taxonomyFile, numberOfEntry)

    for item in self.IdItemMap.values():
      if item.Id != item.ParentId:
        self.ChildToParentMap[item.Id] = item.ParentId
        self.ParentChildTree.setdefault(item.ParentId, []).append(item.Id)

  def GetItem(self, id):
    return(self.IdItemMap[id])

  def GetChildren(self, id):
    result = [self.IdItemMap[childId] for childId in self.ParentChildTree[id]]
    return(result)

  def BuildChildIdParentIdMap(self, rootid):
    result = {}
    subids = set(self.ParentChildTree[rootid])
    stack = deque(self.ParentChildTree[rootid])
    while stack:
      curid = stack.pop()
      if curid in subids:
        parentid = curid
      else:
        parentid = result[curid]

      if curid in self.ParentChildTree:
        curchildrenids = self.ParentChildTree[curid]
        for childid in curchildrenids:
          result[childid] = parentid
          stack.append(childid)
    return(result)

  def BuildChildIdParentIdMapByRank(self, rootid, rank):
    result = {}
    subids = set(id for id in self.ParentChildTree[rootid] if self.IdItemMap[id].Level > level)
    stack = deque(self.ParentChildTree[rootid])
    while stack:
      curid = stack.pop()
      if curid in subids:
        parentid = curid
      else:
        parentid = result[curid]

      if curid in self.ParentChildTree:
        curchildrenids = self.ParentChildTree[curid]
        for childid in curchildrenids:
          result[childid] = parentid
          stack.append(childid)
    return(result)

  def BuildChildIdParentNameMap(self, rootid):
    tmp = self.BuildChildIdParentIdMap(rootid)
    result = {}
    for childid in tmp:
      result[childid] = self.IdItemMap[tmp[childid]].Name
    return(result)

class TestTaxonomyTree(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    super(TestTaxonomyTree, cls).setUpClass()
    cls.tree = TaxonomyTree()
    cls.tree.ReadFromFile("/scratch/cqs_share/references/taxonomy/20200321_taxonomy.txt")

  def testGetChildren(self):
    subitems = self.tree.GetChildren("1")
    self.assertEqual(4, len(subitems))
    self.assertEqual("10239", subitems[0].Id)
    self.assertEqual("12908", subitems[1].Id)
    self.assertEqual("28384", subitems[2].Id)
    self.assertEqual("131567", subitems[3].Id)
    #for subitem in subitems:
    #  print(subitem.Id + "," + subitem.Name)

  def testBuildChildIdParentIdMap(self):
    #get bacteria map
    itemMap = self.tree.BuildChildIdParentIdMap("2")
    self.assertEqual("1224", itemMap["1798413"])

if __name__ == '__main__':
  unittest.main()