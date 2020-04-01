import argparse
import sys
import logging
import os
from datetime import datetime

from .__version__ import __version__
from .database_util import prepare_database, prepare_index
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
  parser_p.add_argument('--maxGenomeInFile', action='store',  type=int, default=500, nargs='?', help='Input number of genome in each output file (default:500)')
  parser_p.add_argument('--prefix', action='store', nargs='?', help='Input prefix of database (default will be date)')
  parser_p.add_argument('-o', '--outputFolder', action='store', nargs='?', help="Output folder", required=NOT_DEBUG)

  # create the parser for the "index" command
  parser_i = subparsers.add_parser('index')
  parser_i.add_argument('-i', '--input', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_i.add_argument('-t', '--thread', action='store', type=int, default=8, nargs='?', help="Thread number")
  parser_i.add_argument('-f', '--force', action='store_true', default=False, help="Ignore existing index file and regenerate all")
  parser_i.add_argument('-s', '--slurm', action='store_true', default=False, help="Generate slurm scripts")
  parser_i.add_argument('-e', '--slurmEmail', action='store', help="Email for slurm status")

  # create the parser for the "bowtie" command
  parser_s = subparsers.add_parser('bowtie')
  parser_s.add_argument('-i', '--input', action='store', nargs='?', help='Input single-end fastq/fasta file', required=NOT_DEBUG)
  parser_s.add_argument('-d', '--databaseListFile', action='store', nargs='?', help='Input database list file', required=NOT_DEBUG)
  parser_s.add_argument('-t', '--thread', action='store', nargs='?', type=int, default=8, help="Thread number")
  parser_s.add_argument('--fastq2fasta', action='store_true', default=False, help="Convert fastq to fasta format for bowtie")
  parser_s.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output file", required=NOT_DEBUG)

  # create the parser for the "count" command
  parser_c = subparsers.add_parser('count')
  parser_c.add_argument('-i', '--input', action='store', nargs='?', help='Input bowtie result list file', required=NOT_DEBUG)
  parser_c.add_argument('-c', '--countFile', action='store', nargs='?', help='Input dupcount list file', required=NOT_DEBUG)
  parser_c.add_argument('-o', '--output', action='store', nargs='?', default="-", help="Output count file", required=NOT_DEBUG)
  
  if not DEBUG and len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()
  #args.command = "database"
  
  if args.command == "database":
    if DEBUG:
      #args.input = "2"
      args.input = "10239"
      args.maxGenomeInFile = 500
      args.outputFolder = "/scratch/cqs_share/references/genbank"
      args.refseq = False

    if args.prefix == None:
      now = datetime.now()
      args.prefix = now.strftime("%Y%m%d_")

    database = "refseq" if args.refseq else "genbank"
    logger = initialize_logger(os.path.join(args.outputFolder, args.prefix + "_spcount_database.log"), args)
    print(args)
    prepare_database(logger, args.input, args.outputFolder, args.maxGenomeInFile, args.prefix, database)
  elif args.command == "index":
    if DEBUG:
      args.input = "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.list"
      args.thread = 8
    logger = initialize_logger(args.input + "_spcount_index.log", args)
    print(args)
    prepare_index(logger, args.input, args.thread, args.force, args.slurm, args.slurmEmail)
  elif args.command == "bowtie":
    if DEBUG:
      args.input = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.gz"
      args.databaseListFile = "/scratch/cqs_share/references/refseq/bacteria/20200321_assembly_summary.txt.files.list"
      args.thread = 32
      args.output = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ.txt"
    print(args)
    if args.fastq2fasta:
      logger = initialize_logger(args.outputPrefix + "_spcount_bowtie_fastq2fasta.log", args)
      bowtie_fastq2fasta(logger, args.input, args.outputPrefix, args.databaseListFile, args.thread)
    else:
      logger = initialize_logger(args.output + "_spcount_bowtie.log", args)
      bowtie(logger, args.input, args.output, args.databaseListFile, args.thread)
  elif args.command == "count":
    if DEBUG:
      args.input = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_bacteria.txt"
      args.countFile = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_clipped_identical.unmapped.fastq.dupcount"
      args.output = "/scratch/cqs/kasey_vickers_projects/testdata/VLDL_WZ_bacteria.count"
    print(args)
    logger = initialize_logger(args.output + "_spcount_count.log", args)
    count(logger, args.input, args.output, args.countFile)
  
if __name__ == "__main__":
  main()
