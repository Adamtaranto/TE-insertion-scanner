#!/usr/bin/env python
#python 2.7.5
#insert-finder.py
#Version 1. Adam Taranto, April 2017
#Contact, Adam Taranto, adam.taranto@anu.edu.au

#########################################################################
# Parse whole genome alignments for signatures of transposon insertion. #                                                        #
#########################################################################

import sys
import os
import argparse
from collections import namedtuple

def dochecks(args):
	# Make outDir if does not exist else set to current dir
	if args.outDir:
		tempPathCheck(args)
		outDir = args.outDir
	else:
		outDir = os.getcwd() 
	return os.path.join(outDir,args.gffOut)

def tempPathCheck(args):
	absOutDir = os.path.abspath(args.outDir)
	if not os.path.isdir(absOutDir):
		os.makedirs(absOutDir)


def readLASTZ(infile,minID=90):
	'''Read in LASTZ result file from LASTZ_genome_align.sh
	Populate nested dictionary of hits keyed by Target and Query scaffold names.'''
	with open(infile) as f:
		content = f.readlines()
	content = [x.strip().split() for x in content] 
	hitsDict = dict()
	hitTup = namedtuple('Elem', ['t_start','t_end','t_strand','q_start','q_end','q_strand','idPct','UID'])
	counter = 0
	# Read in split rows, ignore lines begining with '#'
	for row in content:
		if row[0][0] == "#":
			continue
		elif float(row[9]) >= minID:
			counter += 1
			UID = counter
			t_name	= str(row[0])
			t_strand	= str(row[1])
			t_start	= int(row[2]) - 1
			t_end	= int(row[3]) - 1
			q_name	= str(row[4])
			q_strand	= str(row[5])
			idPct	= float(row[9])
			if int(row[6]) < int(row[7]): 
				q_start	= int(row[6]) - 1
				q_end	= int(row[7]) - 1
			else:
				print("Inverting query sequence coordinates for record number: ", UID)
				q_end	= int(row[6]) - 1
				q_start	= int(row[7]) - 1
			# Create target scaffold dict if not seen
			if t_name not in hitsDict.keys():
				hitsDict[t_name] = dict()
			# Create query scaffold list if not seen
			if q_name not in hitsDict[t_name].keys():
				hitsDict[t_name][q_name] = list()
			# Write record to target:query list as tuple
			hitsDict[t_name][q_name].append(hitTup(t_start,t_end,t_strand,q_start,q_end,q_strand,idPct,UID))
	return hitsDict

def getHitPairs(hits,args):
		valid_elem = []
		for t_name in hits.keys():
				for q_name in hits[t_name].keys():
						for hit in hits[t_name][q_name]:
								for mate in hits[t_name][q_name]:
										if (mate.UID != hit.UID and \
												mate.q_strand == hit.q_strand and \
												abs(mate.idPct - hit.idPct) <= args.maxIdentDiff and \
												mate.t_start - hit.t_end >= args.minInsert and \
												mate.t_start - hit.t_end <= args.maxInsert):
												if (hit.q_strand == '+' and \
														mate.q_start - hit.q_end <= args.qGap and \
														mate.q_start - hit.q_end >= 0-args.maxTSD):
														valid_elem.append(((t_name,q_name),hit,mate))
												elif (hit.q_strand == '-' and \
														hit.q_start - mate.q_end <= args.qGap and \
														hit.q_start - mate.q_end >= 0-args.maxTSD):
														valid_elem.append(((t_name,q_name),hit,mate))
		return valid_elem

def formatGFFline(pair,featureID):
	seqid = str(pair[0][0])
	source = 'InsertFinder'
	feature_type = "Insertion"
	start = pair[1].t_end + 1
	end = pair[2].t_start - 1
	score = '.'
	strand = pair[1].t_strand
	phase = '.'
	leftflank = pair[0][1] + '_' + pair[1].q_strand + '_' + str(pair[1].q_start) + '_' + str(pair[1].q_end)
	rightflank = pair[0][1] + '_' + pair[2].q_strand + '_' + str(pair[2].q_start) + '_' + str(pair[2].q_end)
	attributes  = 'ID=' + str(featureID) + ';len=' + str(end - start) + \
								';leftflank=' + leftflank + ';rightflank=' + rightflank + \
								';leftID=' + str(pair[1].idPct) + ';rightID=' + str(pair[2].idPct) + '\n'
	return '\t'.join([seqid,source,feature_type,str(start),str(end),score,strand,phase,attributes])

def getTSD(pair,parentID):
	TSDlen = None
	if (pair[1].q_strand == '+' and pair[1].q_end > pair[2].q_start):
		TSDlen = pair[1].q_end - pair[2].q_start
	elif (pair[1].q_strand == '-' and pair[1].q_start < pair[2].q_end ):
		TSDlen = pair[2].q_end - pair[1].q_start
	if TSDlen:
		TSDl_start = pair[2].t_start
		TSDl_end = pair[2].t_start + TSDlen
		TSDr_start = pair[1].t_end - TSDlen
		TSDr_end = pair[1].t_end
		strand = pair[1].q_strand
		seqid = str(pair[0][0])
		source = 'InsertFinder'
		feature_type = "TSD"
		attributesR  = 'ID=' + str(parentID) + '_TSD_R' + ';Parent=' + parentID + ';len=' + str(TSDlen) + '\n'
		attributesL  = 'ID=' + str(parentID) + '_TSD_L' + ';Parent=' + parentID + ';len=' + str(TSDlen) + '\n'
		# Yield left
		yield '\t'.join([seqid,source,feature_type,str(TSDl_start),str(TSDl_end),'.',strand,'.',attributesL])
		# Yield right
		yield '\t'.join([seqid,source,feature_type,str(TSDr_start),str(TSDr_end),'.',strand,'.',attributesR])
	else:
		return None

def writeGFFlines(validPairs):
	yield '#gff-version 3\n#seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes\n'
	counter = 0
	fillLen = len(str(abs(len(validPairs)))) + 1
	for pair in validPairs:
		counter += 1
		featureID = 'IS_' + str(counter).zfill(fillLen)
		yield formatGFFline(pair,featureID)
		TSDline = getTSD(pair,featureID)
		if TSDline:
			for x in TSDline:
				yield x

def mainArgs():
		parser = argparse.ArgumentParser(
				description='Parse whole genome alignments for signatures of transposon insertion.',
				prog='insert-finder')
		parser.add_argument('--version',
												action='version',
												version='insert-finder-0.0.1')
		parser.add_argument('--maxTSD',
												type=int,
												default=100,
												help='Maximum overlap of insertion flanking sequences in QUERY genome to be considered as target site duplication.\n \
												Flank pairs with greater overlaps will be discarded \n \
												Note: Setting this value too high may result in tandem duplications in the target genome being falsely classified as insertion events.')
		parser.add_argument('--maxInsert',
												type=int,
												default=20000,
												help='Maximum length of sequence to consider as an insertion event.')
		parser.add_argument('--minInsert',
												type=int,
												default=100,
												help='Minimum length of sequence to consider as an insertion event. \n \
												Note: If too short may detect small non-TE indels.')
		parser.add_argument('--qGap',
												type=int,
												default=100,
												help='Maximum gap allowed between aligned flanks in QUERY sequence. \n \
												Equivalent to target sequence deleted upon insertions event.')
		parser.add_argument('--minIdent',
												type=int,
												default=90,
												help='Minimum identity for a hit to be considered.')
		parser.add_argument('--maxIdentDiff',
												type=float,
												default=20,
												help='Maximum divergence in identity (to query) allowed between insert flanking sequences.')
		parser.add_argument('-i',
												'--infile',
												type=str,
												required=True,
												help='Input file containing tab delimited LASTZ alignment data.')
		parser.add_argument('--outDir',
												default=None,
												help='Optional: Directory to write output to.')
		parser.add_argument('--gffOut',
												type=str,
												default="candidate_insertion_events.gff3",
												help='Write features to this file as gff3')
		args = parser.parse_args()
		return args

def main(args):
	# Check existence of output directory
	outPath = dochecks(args)
	# Read in LASTZ hits file
	hits = readLASTZ(args.infile,minID=args.minIdent)
	# Screen for candidate insertion events
	validPairs 	=  getHitPairs(hits,args)
	# Write insertions and TSDs to gff file
	with open(outPath, 'w') as f:
		for x in writeGFFlines(validPairs):
			f.write(x)

if __name__== '__main__':
	args = mainArgs()
	main(args)