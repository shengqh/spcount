import argparse
import sys
import logging
import os
from datetime import datetime

from .__version__ import __version__
from .taxonomy_util import prepare_taxonomy
from .database_util import download_assembly_summary, prepare_segment_database, prepare_index, fastq_to_database
from .bowtie_util import bowtie, bowtie_fastq2fasta
from .count_util import bowtie_count, count_table
from .visualization_util import krona

def initialize_logger(logfile, args):
  logger = logging.getLogger('spcount')
  loglevel = logging.INFO
  logger.setLevel(loglevel)

  formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)-8s - %(message)s')    
 
  # create console handler and set level to info
  handler = logging.StreamHandler()
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  # create error file handler and set level to error
  handler = logging.FileHandler(logfile, "w")
  handler.setLevel(loglevel)
  handler.setFormatter(formatter)
  logger.addHandler(handler)
 
  return(logger)

def runCommand(command, logger):
  logger.info("run : " + command )
  os.system(command)

def main():
  parser = argparse.ArgumentParser(description="spcount " + __version__,
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  
  DEBUG = False
  NOT_DEBUG = not DEBUG

  now = datetime.now()
  date_time = now.strftime("%Y%m%d_")
  
  subparsers = parser.add_subparsers(dest="command")

  parser_t = subparsers.add_parser('dl_taxonomy')
  parser_t.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)

  parser_s = subparsers.add_parser('dl_assembly_summary')
  parser_s.add_argument('-d', '--database', action='store', default="genbank", help='Input database (genbank or refseq, default is genbank database)')
  parser_s.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)

  parser_segment = subparsers.add_parser('dl_genome')
  parser_segment.add_argument('-i', '--taxonomy_id', action='store', type=int, default=2, nargs='?', required=NOT_DEBUG, help='Input taxonomy id (for example, 2 for bacteria)')
  parser_segment.add_argument('-t', '--taxonomy_file', action='store', nargs='?', required=NOT_DEBUG, help='Input taxonomy file')
  parser_segment.add_argument('-a', '--assembly_summary_file', action='store', nargs='?', required=NOT_DEBUG, help='Input assembly summary file')
  parser_segment.add_argument('-n', '--maximum_genome_in_file', action='store', type=int, default=500, nargs='?', required=NOT_DEBUG, help='Input number of genome in each output file (default:500)')
  parser_segment.add_argument('-p', '--prefix', action='store', nargs='?', required=NOT_DEBUG, help='Input prefix of database')
  parser_segment.add_argument('-r', '--reference_representative_only', action='store_true', required=NOT_DEBUG, help='Use reference or representative genome only')
  parser_segment.add_argument('-o', '--output_folder', action='store', nargs='?', required=NOT_DEBUG, help="Output folder")

  # create the parser for the "index" command
  parser_index = subparsers.add_parser('bowtie_index')
  parser_index.add_argument('-i', '--input', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_index.add_argument('-t', '--thread', action='store', type=int, default=8, nargs='?', help="Thread number")
  parser_index.add_argument('-f', '--force', action='store_true', default=False, help="Ignore existing index file and regenerate all")
  parser_index.add_argument('-s', '--slurm_template', action='store', default=False, help="Input slurm template")

  parser_count = subparsers.add_parser('bowtie_count')
  parser_count.add_argument('-i', '--input', action='store', nargs='?', help='Input BAM list file', required=NOT_DEBUG)
  parser_count.add_argument('-c', '--count', action='store', nargs='?', help='Input count file', required=NOT_DEBUG)
  parser_count.add_argument('-s', '--species', action='store', nargs='?', help='Input species file', required=NOT_DEBUG)
  parser_count.add_argument('-t', '--species_column', action='store', default="species", nargs='?', help='Input species column')
  parser_count.add_argument('-o', '--output', action='store', nargs='?', help="Output summary file", required=NOT_DEBUG)

  parser_table = subparsers.add_parser('count_table')
  parser_table.add_argument('-i', '--input', action='store', nargs='?', help='Input count list file', required=NOT_DEBUG)
  parser_table.add_argument('-t', '--taxonomy_file', action='store', nargs='?', required=NOT_DEBUG, help='Input taxonomy file')
  parser_table.add_argument('-s', '--species', action='store', nargs='?', help='Input species file', required=NOT_DEBUG)
  parser_table.add_argument('-c', '--species_column', action='store', default="species", nargs='?', help='Input species column')
  parser_table.add_argument('-a', '--aggregate_rate', action='store', type=float, default=0.95, help='Input aggregate rate (default 0.95)')
  parser_table.add_argument('-o', '--output_prefix', action='store', nargs='?', help="Output prefix", required=NOT_DEBUG)
  parser_table.add_argument('-d', '--debug_mode', action='store_true', help="Debug mode")

  parser_krona = subparsers.add_parser('krona')
  parser_krona.add_argument('-i', '--input', action='store', nargs='?', help='Input tree count file', required=NOT_DEBUG)
  parser_krona.add_argument('-g', '--group_file', action='store', nargs='?', required=NOT_DEBUG, help='Input group file')
  parser_krona.add_argument('-t', '--taxonomy_folder', action='store', nargs='?', required=NOT_DEBUG, help='Path to directory containing a Krona taxonomy database (taxonomy.tab) to use.')
  parser_krona.add_argument('-o', '--output_prefix', action='store', nargs='?', help="Output prefix", required=NOT_DEBUG)

  # parser_fastq_to_database = subparsers.add_parser('fastq_to_database')
  # parser_fastq_to_database.add_argument('-i', '--input', action='store', nargs='?', help='Input FASTQ file', required=NOT_DEBUG)
  # parser_fastq_to_database.add_argument('--sample_name', action='store', nargs='?', help='Input sample name', required=NOT_DEBUG)
  # parser_fastq_to_database.add_argument('--reads_per_file', action='store',  type=int, default=5000000, nargs='?', help='Input number of reads in each database file')
  # parser_fastq_to_database.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)

  # # create the parser for the "bowtie" command
  # parser_bowtie = subparsers.add_parser('bowtie')
  # parser_bowtie.add_argument('-i', '--input', action='store', nargs='?', help='Input single-end fastq/fasta file', required=NOT_DEBUG)
  # parser_bowtie.add_argument('-d', '--databaseListFile', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  # parser_bowtie.add_argument('-t', '--thread', action='store', nargs='?', type=int, default=8, help="Thread number")
  # parser_bowtie.add_argument('--fastq2fasta', action='store_true', default=False, help="Convert fastq to fasta format for bowtie")
  # parser_bowtie.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file", required=NOT_DEBUG)

  # # create the parser for the "count" command
  # parser_count = subparsers.add_parser('count')
  # parser_count.add_argument('-i', '--input', action='store', nargs='?', help='Input bowtie result list file', required=NOT_DEBUG)
  # parser_count.add_argument('-c', '--countFile', action='store', nargs='?', help='Input dupcount list file', required=NOT_DEBUG)
  # parser_count.add_argument('--category_name', action='store', nargs='?', help="Input category name and ignore the one in bowtie result file")
  # parser_count.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output count file", required=NOT_DEBUG)
  
  # parser_sequential_count = subparsers.add_parser('sequential_count')
  # parser_sequential_count.add_argument('-i', '--input', action='store', nargs='?', help='Input fastq list file', required=NOT_DEBUG)
  # parser_sequential_count.add_argument('-d', '--dbFile', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  # parser_sequential_count.add_argument('-c', '--countFile', action='store', nargs='?', help='Input dupcount list file', required=NOT_DEBUG)
  # parser_sequential_count.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)
  
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  #args.command = "database"
  
  if args.command == "dl_taxonomy":
    logger = initialize_logger(args.output + ".log", args)
    print(args)
    prepare_taxonomy(logger, args.output)
  elif args.command == "dl_assembly_summary":
    logger = initialize_logger(args.output + ".log", args)
    print(args)
    download_assembly_summary(logger, args.output, args.database)
  elif args.command == "dl_genome":
    logger = initialize_logger(os.path.join(args.output_folder, args.prefix + ".log"), args)
    print(args)
    prepare_segment_database(logger, args.taxonomy_file, args.assembly_summary_file, args.taxonomy_id, args.output_folder, args.prefix, args.maximum_genome_in_file, args.reference_representative_only)
  elif args.command == "bowtie_index":
    if DEBUG:
      args.input = "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.list"
      args.thread = 8
    logger = initialize_logger(args.input + ".log", args)
    print(args)
    prepare_index(logger, args.input, args.thread, args.force, args.slurm_template)
  elif args.command == 'bowtie_count':
    logger = initialize_logger(args.output + ".log", args)
    print(args)
    bowtie_count(logger, 
                 input_list_file = args.input, 
                 output_file = args.output,
                 count_file = args.count, 
                 species_file = args.species, 
                 species_column = args.species_column)
  elif args.command == 'count_table':
    logger = initialize_logger(args.output_prefix + ".log", args)
    print(args)
    count_table(logger, 
                input_list_file = args.input, 
                output_prefix = args.output_prefix,
                taxonomy_file = args.taxonomy_file,
                species_file = args.species, 
                species_column = args.species_column,
                aggregate_rate = args.aggregate_rate,
                debug_mode = args.debug_mode)
  elif args.command == "krona":
    logger = initialize_logger(args.output_prefix + ".log", args)
    print(args)
    krona(logger, 
      treeFile = args.input, 
      groupFile = args.group_file, 
      taxonomyFolder = args.tax,
      outputPrefix = args.output_prefix)
  elif args.command == "fastq_to_database":
    if DEBUG:
      #args.input = "2"
      args.input = "/data/vickers_lab/20200828_5059_ES/NextFlex/5059-ES-1_S01_L005_R1_001.fastq.gz"
      args.output = "/scratch/vickers_lab/projects/20200831_5059_ES_scRNA_map_WGS_byTiger/database/5059_ES_1_R1"
      args.sample_name = "5059_ES_1_R1"
      args.reads_per_file = 5000000

    logger = initialize_logger(args.output + ".log", args)
    print(args)
    fastq_to_database(logger, args.input, args.sample_name, args.output, args.reads_per_file)
  elif args.command == "bowtie":
    if DEBUG:
      args.input = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.gz"
      args.databaseListFile = "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.list"
      args.thread = 32
      args.output = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ.txt"
    print(args)
    if args.fastq2fasta:
      logger = initialize_logger(args.outputPrefix + ".log", args)
      bowtie_fastq2fasta(logger, args.input, args.outputPrefix, args.databaseListFile, args.thread)
    else:
      logger = initialize_logger(args.output + ".spcount.log", args)
      bowtie(logger, args.input, args.output, args.databaseListFile, args.thread)
  elif args.command == "count":
    if DEBUG:
      args.input = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_bacteria.txt"
      args.countFile = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.dupcount"
      args.output = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_bacteria.count"
    print(args)
    logger = initialize_logger(args.output + ".log", args)
    count(logger, args.input, args.output, args.countFile, args.category_name)
  elif args.command == "sequential_count":
    logger = initialize_logger(args.output + ".log", args)
    #sequential_count(logger, args.input, args.dbFile, args.output, args.countFile)
  
if __name__ == "__main__":
  main()
