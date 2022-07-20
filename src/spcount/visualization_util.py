import pandas as pd
import subprocess
import logging
import os

def read_taxonomy_name_map(taxonomy_file):
  name_map = {}
  with open(taxonomy_file, "rt") as fin:
    for line in fin:
      parts = line.rstrip().split('\t')
      name_map[parts[4]] = {
        "id":int(parts[0]),
        "rank_level": int(parts[1]),
        "parent_id": int(parts[2]),
        "rank": parts[3],
        "name": parts[4]
      }
  return(name_map)

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

def get_taxonomy_id(name_map, row):  
  return (name_map[row['Feature']]['id'])

def get_rank(name_map, row):  
  return (name_map[row['Feature']]['rank'])

def krona(logger, treeFile, groupFile, taxonomyFolder, outputPrefix):
  tree_data = pd.read_csv(treeFile, sep="\t")
  if tree_data.columns[1] != "TaxonomyId":
    name_map = read_taxonomy_name_map(os.path.join(taxonomyFolder, "taxonomy.tab"))
    tree_data.insert(loc=1, column="TaxonomyId", value=tree_data.apply(lambda row: get_taxonomy_id(name_map, row), axis=1))
    tree_data.insert(loc=2, column="Rank", value=tree_data.apply(lambda row: get_rank(name_map, row), axis=1))
    treeFile=outputPrefix + ".tree.count"
    tree_data.to_csv(treeFile, sep="\t", index=False)

  logger.info("Start sample krona ...")
  
  draw_krona(logger, treeFile, taxonomyFolder, outputPrefix)

  logger.info("Start group krona ...")

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

  group_file = outputPrefix + ".group.txt"
  group_data.to_csv(group_file, sep="\t", index=None)

  draw_krona(logger, group_file, taxonomyFolder, outputPrefix)

  logger.info("done")

if __name__ == "__main__":
  logger = logging.getLogger('spcount')
  krona(logger, "/scratch/cqs/shengq2/spcount/RA_97_93.tree.count", "/scratch/cqs/shengq2/spcount/fileList2.txt", "/data/cqs/references/spcount/", "/scratch/cqs/shengq2/spcount/RA_97_93")
