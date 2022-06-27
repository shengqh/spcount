import pandas as pd
import numpy as np
import logging

def shannon_value(x):
    return -(x * np.log(x))
  
def shannon(logger, countFile, outputFile):
  logger.info("read " + countFile + "...")
  count_data = pd.read_csv(countFile, sep="\t", index_col=0)

  logger.info("calculate shannon diversity index ...")
  perc_count = count_data.div(count_data.sum(axis=0), axis=1)

  s_values = perc_count.applymap(shannon_value)
  shannon_values=s_values.sum(axis=0)
  shannon_values.name = "shannon"
  shannon_values.index.name = "sample"

  logger.info(f"save to {outputFile} ...")
  shannon_values.to_csv(outputFile)
  
  logger.info("done")

if __name__ == "__main__":
  logger = logging.getLogger('spcount')
  shannon(logger, "/scratch/vickers_lab/projects/20220509_smRNA_6254_AC_samples_1_13_redo/refseq_bacteria_table/smRNA_6254_AC_redo.species.estimated.count", "/scratch/vickers_lab/projects/20220509_smRNA_6254_AC_samples_1_13_redo/refseq_bacteria_table/smRNA_6254_AC_redo.shannon.csv")
