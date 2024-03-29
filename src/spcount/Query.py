import pandas as pd

class Query(object):
  __slots__ = "sample", "name", "count", "species_list", "rank", "estimated_count", 

  def __init__(self, sample, name, count, species_list, rank = "species"):
    self.sample = sample
    self.name = name
    self.count = count
    self.species_list = species_list
    self.rank = rank
    self.estimated_count = 0
  
  def remove_species(self, species_list):
    for species in species_list:
      if species in self.species_list:
          self.species_list.remove(species)

  def estimate_count(self):
    self.estimated_count = self.count / len(self.species_list)

  def get_unique_rank(self, species_taxonomy_map, ranks=['genus', 'family', 'order', 'class', 'phylum', 'superkingdom'], aggregate_rate=0.95):
    if len(self.species_list) == 1:
      return self.rank, self.species_list[0]
    else:
      for rank in ranks:
        rank_list = [species_taxonomy_map[s][rank] for s in self.species_list]
        vc = pd.Series(rank_list).value_counts(sort=True, ascending=False)
        ar = vc[0] / len(rank_list)
        if ar >= aggregate_rate:
          rank_name = vc.keys()[0]
          if rank_name == "Unclassified":
            continue
          return rank, rank_name
      raise Exception(f"Cannot find aggregated rank for {'.'.join(self.species_list)}")
  
  def aggregate_to_rank(self, species_taxonomy_map, rank, aggregate_rate=0.95):
    if self.rank == rank:
      if len(self.species_list) > 1:
        return("AmbiguousRanks")
      else:
        return(self.species_list[0])

    rank_list = [species_taxonomy_map[s][rank] for s in self.species_list]
    vc = pd.Series(rank_list).value_counts(sort=True, ascending=False)
    ar = vc.iloc[0] / len(rank_list)
    if ar >= aggregate_rate:
      rank_name = vc.keys()[0]
      return(rank_name)
    else:
      return("AmbiguousRanks")
  
