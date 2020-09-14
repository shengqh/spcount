import argparse
import sys
import logging
import os
from datetime import datetime

from .__version__ import __version__
from .database_util import prepare_database, prepare_database_as_whole, prepare_index, fastq_to_database
from .bowtie_util import bowtie, bowtie_fastq2fasta
from .count_util import count

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

  # create the parser for the "database" command
  parser_p = subparsers.add_parser('database')
  parser_p.add_argument('-i', '--input', action='store', nargs='?', help='Input root taxonomy id (2 for bacteria)', required=NOT_DEBUG)
  parser_p.add_argument('--refseq', action='store_true', help='Use refseq database (default is genbank database)')
  parser_p.add_argument('--maximum_genome_in_file', action='store',  type=int, default=500, nargs='?', help='Input number of genome in each output file (default:500)')
  parser_p.add_argument('--prefix', action='store', nargs='?', help='Input prefix of database (default will be date)')
  parser_p.add_argument('-o', '--output_folder', action='store', nargs='?', help="Output folder", required=NOT_DEBUG)

  # create the parser for the "database_whole" command
  parser_w = subparsers.add_parser('database_whole')
  parser_w.add_argument('-i', '--input', action='store', nargs='?', help='Input root taxonomy id (2 for bacteria)', required=NOT_DEBUG)
  parser_w.add_argument('--refseq', action='store_true', help='Use refseq database (default is genbank database)')
  parser_w.add_argument('--prefix', action='store', nargs='?', help='Input prefix of database (default will be date)')
  parser_w.add_argument('-o', '--output_folder', action='store', nargs='?', help="Output folder", required=NOT_DEBUG)

  parser_fastq_to_database = subparsers.add_parser('fastq_to_database')
  parser_fastq_to_database.add_argument('-i', '--input', action='store', nargs='?', help='Input FASTQ file', required=NOT_DEBUG)
  parser_fastq_to_database.add_argument('--sample_name', action='store', nargs='?', help='Input sample name', required=NOT_DEBUG)
  parser_fastq_to_database.add_argument('--reads_per_file', action='store',  type=int, default=5000000, nargs='?', help='Input number of reads in each database file')
  parser_fastq_to_database.add_argument('-o', '--output', action='store', nargs='?', help="Output file", required=NOT_DEBUG)

  # create the parser for the "index" command
  parser_index = subparsers.add_parser('index')
  parser_index.add_argument('-i', '--input', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_index.add_argument('-t', '--thread', action='store', type=int, default=8, nargs='?', help="Thread number")
  parser_index.add_argument('-f', '--force', action='store_true', default=False, help="Ignore existing index file and regenerate all")
  parser_index.add_argument('-s', '--slurm', action='store_true', default=False, help="Generate slurm scripts")
  parser_index.add_argument('-e', '--slurm_email', action='store', help="Email for slurm status")

  # create the parser for the "bowtie" command
  parser_bowtie = subparsers.add_parser('bowtie')
  parser_bowtie.add_argument('-i', '--input', action='store', nargs='?', help='Input single-end fastq/fasta file', required=NOT_DEBUG)
  parser_bowtie.add_argument('-d', '--databaseListFile', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_bowtie.add_argument('-t', '--thread', action='store', nargs='?', type=int, default=8, help="Thread number")
  parser_bowtie.add_argument('--fastq2fasta', action='store_true', default=False, help="Convert fastq to fasta format for bowtie")
  parser_bowtie.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file", required=NOT_DEBUG)

  # create the parser for the "count" command
  parser_count = subparsers.add_parser('count')
  parser_count.add_argument('-i', '--input', action='store', nargs='?', help='Input bowtie result list file', required=NOT_DEBUG)
  parser_count.add_argument('-c', '--countFile', action='store', nargs='?', help='Input dupcount list file', required=NOT_DEBUG)
  parser_count.add_argument('--category_name', action='store', nargs='?', help="Input category name and ignore the one in bowtie result file")
  parser_count.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output count file", required=NOT_DEBUG)
  
  parser_sequential_count = subparsers.add_parser('sequential_count')
  parser_sequential_count.add_argument('-i', '--input', action='store', nargs='?', help='Input fastq list file', required=NOT_DEBUG)
  parser_sequential_count.add_argument('-d', '--dbFile', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_sequential_count.add_argument('-c', '--countFile', action='store', nargs='?', help='Input dupcount list file', required=NOT_DEBUG)
  parser_sequential_count.add_argument('-o', '--output', action='store', nargs='?', help="Output count file", required=NOT_DEBUG)
  
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  #args.command = "database"
  
  if args.command == "database":
    if DEBUG:
      #args.input = "2"
      args.input = "10239"
      args.maximum_genome_in_file = 500
      args.output_folder = "/scratch/cqs_share/references/genbank"
      args.refseq = False

    if args.prefix == None:
      now = datetime.now()
      args.prefix = now.strftime("%Y%m%d_")

    database = "refseq" if args.refseq else "genbank"
    logger = initialize_logger(os.path.join(args.output_folder, args.prefix + ".log"), args)
    print(args)
    prepare_database(logger, args.input, args.output_folder, args.maximum_genome_in_file, args.prefix, database)
  elif args.command == "database_whole":
    if DEBUG:
      #args.input = "2"
      args.input = "10239"
      args.output_folder = "/scratch/cqs_share/references/genbank"
      args.refseq = False

    if args.prefix == None:
      now = datetime.now()
      args.prefix = now.strftime("%Y%m%d_")

    database = "refseq" if args.refseq else "genbank"
    logger = initialize_logger(os.path.join(args.output_folder, args.prefix + ".log"), args)
    print(args)
    prepare_database_as_whole(logger, args.input, args.output_folder, args.prefix, database)
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
  elif args.command == "index":
    if DEBUG:
      args.input = "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.list"
      args.thread = 8
    logger = initialize_logger(args.input + ".log", args)
    print(args)
    prepare_index(logger, args.input, args.thread, args.force, args.slurm, args.slurm_email)
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
    sequential_count(logger, args.input, args.dbFile, args.output, args.countFile)

  
if __name__ == "__main__":
  main()
