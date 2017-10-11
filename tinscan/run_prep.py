#!/usr/bin/env python
import os 
import sys
import argparse
from Bio import SeqIO
from shlex import quote

def mainArgs():
	parser = argparse.ArgumentParser(description='Split multifasta genome files into directories for A and B genomes.',prog='tinscan-prep')
	# Input options
	parser.add_argument('-A','--target',type=str,required=True,help='Multifasta containing A genome.')
	parser.add_argument('-B','--query',type=str,required=True,help='Multifasta containing B genome.')
	# Output options
	parser.add_argument('--adir',type=str,default=None,help='A genome sub-directory within outdir')
	parser.add_argument('--bdir',type=str,default=None,help='B genome sub-directory within outdir')
	parser.add_argument('-d','--outdir',type=str,default=None,help='Write split directories within this directory. (Default: cwd)')
	args = parser.parse_args()
	return args

def	splitFasta(infile,outdir,unique=True):
	seen = list()
	for rec in SeqIO.parse(infile, "fasta"):
		if str(rec.id) in seen and unique:
			print("Non-unique name in genome: %s. Quitting." % str(rec.id))
			sys.exit(1)
		else:
			seen.append(str(rec.id))
		outfile	= os.path.join(outdir, rec.id + ".fa")
		with open(outfile, "w") as handle:
			SeqIO.write(rec, handle, "fasta")

def isfile(path):
	path = os.path.abspath(path)
	if not os.path.isfile(path):
		print("Input file not found: %s" % quote(path))
		sys.exit(1)
	else:
		return path

def check_paths(args):
	# Check for genome files
	A_genome = isfile(args.target)
	B_genome = isfile(args.query)
	# Check and set outdir
	if args.outdir:
		outdir = os.path.abspath(args.outdir)
		if not os.path.isdir(args.outdir):
			os.makedirs(outdir)
	else:
		outdir = os.getcwd() 
	# Set split directory paths
	if not args.adir:
		A_dir = os.path.join(outdir,"A_target_split")
	else:
		A_dir = os.path.join(outdir,args.adir)
	if not args.bdir:
		B_dir = os.path.join(outdir,"B_query_split")
	else:
		B_dir = os.path.join(outdir,args.bdir)
	# Create dirs if do not exist
	if not os.path.isdir(A_dir):
			os.makedirs(A_dir)
	if not os.path.isdir(B_dir):
			os.makedirs(B_dir)
	# Return paths
	return A_genome,B_genome,A_dir,B_dir

def main():
	# Get cmd line args
	args = mainArgs()
	# Check files and set outpaths
	A_genome,B_genome,A_dir,B_dir = check_paths(args)
	# Read and split input genomes - enforce unique names within genomes
	splitFasta(A_genome,A_dir)
	splitFasta(B_genome,B_dir)