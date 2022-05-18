import logging
import gzip

def take_count(elem):
    return elem[0]

def collapse_fastq(logger, inputFile, outputFilePrefix):
  qname_map = {}

  logger.info(f"reading {inputFile} ...")
  fin = gzip.open(inputFile, "rt") if inputFile.endswith(".gz") else open(inputFile, "rt")
  with fin:
    icount = 0
    while True:
      query = fin.readline()
      if not query:
        break
      sequence = fin.readline()
      line2 = fin.readline()
      score = fin.readline()

      if query == '':
        raise Exception(f"Error format: query={query.rstrip()}")

      if len(sequence) != len(score):
        raise Exception(f"Error format: seq={sequence.rstrip()}, score={score.rtrip()}")

      icount += 1
      if icount % 100000 == 0:
        logger.info(icount)
        #break

      seq_array = qname_map.get(sequence)
      if seq_array == None:
        qname = query.split('\t', 1)[0].split(' ', 1)[0]
        seq_array = [1, qname, sequence, line2, score]
        qname_map[sequence] = seq_array
      else:
        seq_array[0] = seq_array[0] + 1

  logger.info(f"writing {outputFilePrefix} ...")
  queries = list(qname_map.values())
  queries.sort(key=take_count, reverse=True)
  with gzip.open(outputFilePrefix + ".gz", "wt") as fout:
    with open(outputFilePrefix + ".dupcount", "wt") as fcount:
      fcount.write("Query\tCount\tSequence\n")
      for query in queries:
        fcount.write(f"{query[1][1:]}\t{query[0]}\t{query[2]}")
        fout.write(f"{query[1]}\n{query[2]}{query[3]}{query[4]}")

  logger.info("done")

if __name__ == "__main__":
  logger = logging.getLogger('collapse_fastq')
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')
  collapse_fastq(logger, 
            "/scratch/vickers_lab/projects/20220502_8017_AC_smRNA_mouse_byMars/intermediate_data/cutadapt/result/LD-C1_clipped.fastq.gz", 
            "/scratch/vickers_lab/projects/20220502_8017_AC_smRNA_mouse_byMars/preprocessing/identical/result/test.fastq")
