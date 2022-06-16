import pandas as pd
import subprocess
import logging

def draw_krona(logger, treeFile, taxonomyFolder, outputPrefix):
  tree_data = pd.read_csv(treeFile, sep="\t")
  for ind in range(3, tree_data.shape[1]):
    sample = tree_data.columns[ind]
    logger.info(f"processing {sample} ...")
    
    count_file = f"{outputPrefix}.{sample}.txt"
    count_df = tree_data.iloc[:,[0,1,ind]]
    count_df = count_df[count_df.iloc[:,2] > 0]
    count_df.to_csv(count_file, sep="\t", index=None, header=None)
    
    args = ['ktImportTaxonomy', '-o', f"{outputPrefix}.{sample}.html", 
                     '-tax', taxonomyFolder, 
                     '-m', '3',
                     f"{count_file},{sample}"]
    logger.info(" ".join(args))
    subprocess.call(args)

def krona(logger, treeFile, groupFile, taxonomyFolder, outputPrefix):
  logger.info("Start sample krona ...")
  
  draw_krona(logger, treeFile, taxonomyFolder, outputPrefix)

  logger.info("Start group krona ...")

  tree_data = pd.read_csv(treeFile, sep="\t")

  logger.info(f"Read group info {groupFile} ...")
  groups_df=pd.read_csv(groupFile, sep="\t", header=None)
  groups_df.rename(columns={groups_df.columns[0]:"sample", groups_df.columns.values[1]:"group"}, inplace=True)

  group_data=tree_data.iloc[:, [0,1,2]]

  group_data['All'] = tree_data.iloc[:,3:].sum(axis=1)
  for gf in groups_df.groupby('group'):
    gname = gf[0]
    gdf = gf[1]
    gsamples = gdf['sample'].tolist()
    group_data[gname] = tree_data[gsamples].sum(axis=1)

  group_file = outputPrefix + ".group.txt";
  group_data.to_csv(group_file, sep="\t", index=None)

  draw_krona(logger, group_file, taxonomyFolder, outputPrefix)

  logger.info("done")

if __name__ == "__main__":
  logger = logging.getLogger('spcount')
  krona(logger, "/scratch/cqs/shengq2/spcount/RA_97_93.tree.count", "/scratch/cqs/shengq2/spcount/fileList2.txt", "/data/cqs/references/spcount/", "/scratch/cqs/shengq2/spcount/RA_97_93")
