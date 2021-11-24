class CategoryEntry(object):
  def __init__(self, name):
    self.Name = name
    self.TotalCount = 0
    self.Items = []
    self.QueryNames = set()
