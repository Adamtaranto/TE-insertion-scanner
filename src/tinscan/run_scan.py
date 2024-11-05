import argparse
import os
import sys

import tinscan as ts


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Parse whole genome alignments for signatures of transposon insertion.",
        prog="tinscan-find",
    )
    # Input
    parser.add_argument(
        "-i",
        "--infile",
        type=str,
        required=True,
        help="Input file containing tab delimited LASTZ alignment data.",
    )
    # Output
    parser.add_argument(
        "--outdir", default=None, help="Optional: Directory to write output to."
    )
    parser.add_argument(
        "--gffOut",
        type=str,
        default="candidate_insertion_events.gff3",
        help="Write features to this file as gff3.",
    )
    parser.add_argument(
        "--noflanks",
        action="store_false",
        default=True,
        help="If set, do not report flanking hit regions in GFF.",
    )
    # Insert scan settings
    parser.add_argument(
        "--maxTSD",
        type=int,
        default=100,
        help="Maximum overlap of insertion flanking sequences in QUERY genome to be considered as target site duplication. Flank pairs with greater overlaps will be discarded Note: Setting this value too high may result in tandem duplications in the target genome being falsely classified as insertion events.",
    )
    parser.add_argument(
        "--maxInsert",
        type=int,
        default=100000,
        help="Maximum length of sequence to consider as an insertion event.",
    )
    parser.add_argument(
        "--minInsert",
        type=int,
        default=100,
        help="Minimum length of sequence to consider as an insertion event. Note: If too short may detect small non-TE indels.",
    )
    parser.add_argument(
        "--qGap",
        type=int,
        default=100,
        help="Maximum gap allowed between aligned flanks in QUERY sequence. Equivalent to target sequence deleted upon insertion event.",
    )
    parser.add_argument(
        "--minIdent",
        type=int,
        default=90,
        help="Minimum identity for a hit to be considered.",
    )
    parser.add_argument(
        "--maxIdentDiff",
        type=float,
        default=20,
        help="Maximum divergence in identity (to query) allowed between insert flanking sequences.",
    )
    args = parser.parse_args()
    return args


def set_paths(args):
    # Check alignment results file exists
    if not os.path.isfile(args.infile):
        print("Alignment results not found: %s" % args.infile)
        sys.exit(1)
    # Set outdir
    if args.outdir:
        outdir = os.path.abspath(args.outdir)
        if not os.path.isdir(args.outdir):
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Compose path to outfile
    gffout = os.path.join(outdir, args.gffOut)
    return gffout


def main():
    # Get args
    args = mainArgs()
    # Check existence of output directory
    gffout = set_paths(args)
    # Read in LASTZ hits file
    hits = ts.readLASTZ(args.infile, minID=args.minIdent)
    # Screen for candidate insertion events
    validPairs = ts.getHitPairs(hits, args)
    # Write insertions and TSDs to gff file
    with open(gffout, "w") as f:
        for x in ts.writeGFFlines(validPairs, args.noflanks):
            f.write(x)
