class BowtieIndexItem(object):
  def __init__(self, index, category, fasta):
    self.Index = index
    self.Category = category
    self.Fasta = fasta

def readBowtieIndexList(fileName):
  result = []
  with open(fileName, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      result.append(BowtieIndexItem(parts[0], parts[1], parts[2]))
  return(result)

def writeBowtieIndexList(fileName, bowtieIndecies):
  with open(fileName, "wt") as fout:
    fout.write("BowtieIndex\tCategory\tFasta\n")
    for bowtieIndex in bowtieIndecies:
      fout.write("%s\t%s\t%s\n" % (bowtieIndex.Index, bowtieIndex.Category, bowtieIndex.Fasta))

