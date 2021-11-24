from context import spcount

from spcount.Taxonomy import TaxonomyItem, TaxonomyTree

import unittest

class BasicTestSuite(unittest.TestCase):
  def test_ReadFromFile(self):
    tt = TaxonomyTree()
    tt.ReadFromFile("/data/cqs/references/spcount/20211119_taxonomy.txt")
    az = tt.GetItem(6)
    self.assertEqual(335928, az.ParentId)
    self.assertEqual(1224, az.RankMap['phylum'])
    self.assertEqual(28211, az.RankMap['class'])

if __name__ == '__main__':
    unittest.main()