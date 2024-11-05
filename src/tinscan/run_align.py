import argparse
import os
import shutil
import sys

import tinscan


def mainArgs():
    parser = argparse.ArgumentParser(
        description="Align B genome (query) sequences onto A genome (target) using LASTZ.",
        prog="tinscan-align",
    )
    # Input options
    parser.add_argument(
        "--adir",
        type=str,
        required=True,
        help="Name of directory containing sequences from A genome.",
    )
    parser.add_argument(
        "--bdir",
        type=str,
        required=True,
        help="Name of directory containing sequences from B genome.",
    )
    parser.add_argument(
        "--pairs",
        type=str,
        default=None,
        help="Optional: Tab-delimited 2-col file specifying target:query sequence pairs to be aligned",
    )
    # Output options
    parser.add_argument(
        "-d",
        "--outdir",
        type=str,
        default=None,
        help="Write output files to this directory. (Default: cwd)",
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default="tinscan_alignment.tab",
        help="Name of alignment result file.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="If set report LASTZ progress.",
    )
    # LASTZ options
    parser.add_argument(
        "--lzpath",
        type=str,
        default="lastz",
        help="Custom path to LASTZ executable if not in $PATH.",
    )
    parser.add_argument(
        "--minIdt", type=int, default=60, help="Minimum alignment identity to report."
    )
    parser.add_argument(
        "--minLen", type=int, default=100, help="Minimum alignment length to report."
    )
    parser.add_argument(
        "--hspthresh",
        type=int,
        default=3000,
        help="LASTZ min HSP threshold. Increase for stricter matches.",
    )
    args = parser.parse_args()
    return args


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def set_paths(args):
    # Check that A/B genome directories exist
    if not os.path.isdir(args.adir):
        print("Target sequence directory not found: %s" % args.adir)
        sys.exit(1)
    else:
        adir = os.path.abspath(args.adir)
    if not os.path.isdir(args.bdir):
        print("Query sequence directory not found: %s" % args.bdir)
        sys.exit(1)
    else:
        bdir = os.path.abspath(args.bdir)
    # Set outdir
    if args.outdir:
        outdir = os.path.abspath(args.outdir)
        if not os.path.isdir(args.outdir):
            os.makedirs(outdir)
    else:
        outdir = os.getcwd()
    # Compose path to outfile
    outtab = os.path.join(outdir, args.outfile)
    return adir, bdir, outdir, outtab


def main():
    """Do the work."""
    # Get cmd line args
    args = mainArgs()
    # Check for LASTZ
    if missing_tool(args.lzpath):
        print("LASTZ executable was not found at: %s \n Quitting." % args.lzpath)
        sys.exit(1)
    # Set output paths
    adir_path, bdir_path, outdir, outtab = set_paths(args)
    # Import file names to be compaired if set
    if args.pairs and os.path.isfile(args.pairs):
        pairs = tinscan.import_pairs(
            file=os.path.abspath(args.pairs), Adir=adir_path, Bdir=bdir_path
        )
    # Else run all pairwise alignments between A and B genomes
    else:
        pairs = tinscan.get_all_pairs(Adir=adir_path, Bdir=bdir_path)
    # Compose alignment commands
    cmds = tinscan.LASTZ_cmds(
        lzpath=args.lzpath,
        pairs=pairs,
        minIdt=args.minIdt,
        minLen=args.minLen,
        hspthresh=args.hspthresh,
        outfile=outtab,
        verbose=args.verbose,
    )
    # Run alignments
    tinscan.run_cmd(cmds, verbose=args.verbose)
    print("Finished!")
