
class Query(object):
  def __init__(self, name, count, species_list):
    self.name = name
    self.count = count
    self.species_list = species_list
  
  def remove_species(self, species_list):
    for species in species_list:
      if species in self.species_list:
          self.species_list.remove(species)

  def estimate_count(self):
    self.estimated_count = self.count / len(self.species_list)
