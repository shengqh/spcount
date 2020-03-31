import unittest
import requests
import os
import shutil
import tarfile
from collections import deque

class TaxonomyItem(object):
  def __init__(self, name, id, parentid="-1", level=-1):
    self.Name = name
    self.Id = id
    self.Parentid = parentid
    self.Level = level

class TaxonomyTree(object):
  def __init__(self):
    self.IdItemMap = {}
    self.ChildToParentMap = {}
    self.ParentChildTree = {}

  def ReadFromFile(self, taxonomyFile):
    self.IdItemMap.clear()
    self.ChildToParentMap.clear()
    self.ParentChildTree.clear()

    with open(taxonomyFile, "rt") as fin:
      fin.readline()
      for line in fin:
        parts = line.rstrip().split('\t')
        id = parts[0]
        parentId = parts[1]
        name = parts[2]
        item = TaxonomyItem(name, id, parentId)
        self.IdItemMap[id] = item
        self.ChildToParentMap[id] = parentId
        if id != parentId:
          self.ParentChildTree.setdefault(parentId, []).append(id)

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