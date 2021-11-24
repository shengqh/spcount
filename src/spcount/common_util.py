def readFileMap(fileName):
  result = {}
  with open(fileName, "rt") as fin:
    #fin.readline()
    for line in fin:
      parts = line.rstrip().split('\t')
      result[parts[1]] = parts[0]
  return(result)
