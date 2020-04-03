class GenomeItem :
  def __init__(self, category, taxonomyId, accession, name, assemblyLevel, urlFile, assemblyLevelIndex=0):
    self.Category = category
    self.TaxonomyId = taxonomyId
    self.Accession = accession
    self.Name = name
    self.AssemblyLevel = assemblyLevel
    self.UrlFile = urlFile
    self.AssemblyLevelIndex = assemblyLevelIndex

def writeGenomeItems(fileName, items):
  with open(fileName, "wt") as fout:
    fout.write("Category\tTaxonomyId\tAccession\tName\tAssemblyLevel\tUrlFile\n")
    for gi in items:
      fout.write("%s\t%d\t%s\t%s\t%s\t%s\n" % ( gi.Category, gi.TaxonomyId, gi.Accession, gi.Name, gi.AssemblyLevel, gi.UrlFile))

def readGenomeItems(fileName):
  result = []
  with open(fileName, "rt") as fin:
    fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      gi = GenomeItem(parts[0], int(parts[1]), parts[2], parts[3], parts[4], parts[5])
      result.append(gi)
  return(result)
  