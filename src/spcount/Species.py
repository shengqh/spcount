
from collections import defaultdict
from .Query import Query

class Species(object):
  __slots__ = "name", "rank", "taxid", "queries", "queries_set", "is_subset", "is_identical", "identical_species", "sample_query_count", "query_count", "sub_species", "sample_estimated_count", "estimated_count"

  def __init__(self, name, rank=""):
    self.name = name
    self.rank = rank
    self.taxid = ""
    self.queries = defaultdict(list)
    self.queries_set = defaultdict(set)
    self.is_subset = False
    self.is_identical = False
    self.identical_species = []
    self.sample_query_count = {}
    self.sub_species = []
  
  def add_query(self, query):
    self.queries[query.sample].append(query)
    self.queries_set[query.sample].add(query.name)

  def sum_query_count(self):
    self.sample_query_count = {sample:sum([q.count for q in self.queries[sample]]) for sample in self.queries.keys()}
    self.query_count = sum(self.sample_query_count.values())

  def get_query_count_str(self, samples):
    result = "\t".join(str(self.sample_query_count[sample]) if sample in self.sample_query_count else "0" for sample in samples)
    return(result)

  def get_estimated_count_str(self, samples):
    result = "\t".join("{:.2f}".format(self.sample_estimated_count[sample]) if sample in self.sample_estimated_count else "0" for sample in samples)
    return(result)

  def num_of_species(self):
    result = len(self.identical_species) + 1
    return(result)

  def sum_estimated_count(self):
    nos = self.num_of_species()
    self.sample_estimated_count = {sample:nos * sum([q.estimated_count for q in self.queries[sample]]) for sample in self.queries.keys()}
    self.estimated_count = sum(self.sample_estimated_count.values())

  def contains(self, another):
    for sample in another.sample_query_count.keys():
      if not sample in self.sample_query_count:
        return(False)
      if self.sample_query_count[sample] < another.sample_query_count[sample]:
        return(False)

    for sample in another.queries.keys():
      another_q = another.queries_set[sample]
      self_q = self.queries_set[sample]
      if not another_q.issubset(self_q):
          return(False)
      
    return(True)
