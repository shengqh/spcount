from context import spcount

from spcount.count_util import BowtieCountItem, getQueryMap, getCategoryMap, getCategoryUniqueCountMap, getCategoryCountMap, assignCount

import unittest

class BasicTestSuite(unittest.TestCase):
  def test_getQueryMap(self):
    item1 = BowtieCountItem("query1", '+', "chromosome1", 1, "category1", "sample")
    item2 = BowtieCountItem("query1", '+', "chromosome2", 1, "category2", "sample")
    item3 = BowtieCountItem("query2", '+', "chromosome3", 1, "category2", "sample")
    actual = getQueryMap([item1, item2, item3])
    expect = {"query1":[item1, item2], "query2":[item3]}
    self.assertDictEqual(expect, actual)

  def test_getCategoryMap(self):
    item1 = BowtieCountItem("query1", '+', "chromosome1", 1, "category1", "sample")
    item2 = BowtieCountItem("query1", '+', "chromosome2", 1, "category2", "sample")
    item3 = BowtieCountItem("query2", '+', "chromosome3", 1, "category2", "sample")
    actual = getCategoryMap([item1, item2, item3])
    expect = {"category1":[item1], "category2":[item2, item3]}
    self.assertDictEqual(expect, actual)

  def test_getCategoryUniqueCountMap(self):
    item1 = BowtieCountItem("query1", '+', "chromosome1", 1, "category1", "sample", 1)
    item2 = BowtieCountItem("query1", '+', "chromosome2", 1, "category2", "sample", 10)
    item3 = BowtieCountItem("query2", '+', "chromosome3", 1, "category2", "sample", 100)
    actual = getCategoryUniqueCountMap([item1, item2, item3])
    expect = {"category1":1, "category2":2}
    self.assertDictEqual(expect, actual)

  def test_getCategoryCountMap(self):
    item1 = BowtieCountItem("query1", '+', "chromosome1", 1, "category1", "sample", 1)
    item2 = BowtieCountItem("query1", '+', "chromosome2", 1, "category2", "sample", 10)
    item3 = BowtieCountItem("query2", '+', "chromosome3", 1, "category2", "sample", 100)
    actual = getCategoryCountMap([item1, item2, item3])
    expect = {"category1":1, "category2":110}
    self.assertDictEqual(expect, actual)

  def test_assignCount(self):
    item1 = BowtieCountItem("query1", '+', "chromosome1", 1, "category1", "sample")
    item2 = BowtieCountItem("query1", '+', "chromosome2", 1, "category2", "sample")
    item3 = BowtieCountItem("query2", '+', "chromosome3", 1, "category2", "sample")
    items = [item1, item2, item3]
    queryMap = getQueryMap(items)
    assignCount(queryMap, items)
    self.assertEqual(1/3, item1.Count)
    self.assertEqual(2/3, item2.Count)
    self.assertEqual(1, item3.Count)
    

if __name__ == '__main__':
    unittest.main()