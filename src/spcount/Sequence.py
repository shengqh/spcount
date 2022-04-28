class Sequence(object):
  def __init__(self, seq):
    self.seq = seq
    self.sample_query_count = {}
    self.query_count = 0

  def add_sample_count(self, sample, count):
    self.sample_query_count[sample] = count
    self.query_count = self.query_count + count

  def get_sample_count_str(self, samples):
    result = "\t".join(str(self.sample_query_count[sample]) if sample in self.sample_query_count else "0" for sample in samples)
    return(result)
