from context import spcount

from spcount.Taxonomy import load_taxonomy_map

import unittest

class BasicTestSuite(unittest.TestCase):
  def test_load_taxonomy_map(self):
    idItemMap = load_taxonomy_map("/data/cqs/references/spcount/20211119_taxonomy.txt", 100)
    az = idItemMap[6]
    self.assertEqual(335928, az.ParentId)
    self.assertEqual(1224, az.RankMap['phylum'])
    self.assertEqual(28211, az.RankMap['class'])

if __name__ == '__main__':
    unittest.main()