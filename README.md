# Tinscan: TE-Insertion-Scanner

Scan whole genome alignments for transposon insertion signatures.


# Table of contents

* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing Tinscan](#installing-tinscan)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)
      *  [tinscan-prep](#tinscan-prep)
      *  [tinscan-align](#tinscan-align)
      *  [tinscan-find](#tinscan-find)


# Algorithm overview

  1. Perform gapped and chained whole genome alignment of query genome B onto target genome A.
  2. Where two aligned segments are contiguous in B (or separated by no more than --qGap), and 
  3. Separated by an insertion in the range --minInsert:--maxInsert in A, and
  4. At least one flanking alignment in A statisfies the threshold --minIdent, and differs from its mate by no more than --maxIdentDiff %
  5. Log flanks and candidate insertion
  6. Attempt to infer TSDs from internal overlap of flanking alignments in B genome.


# Options and usage 

### Installing Tinscan

Requirements: 
  * [LASTZ](http://www.bx.psu.edu/~rsharris/lastz/) genome alignment tool from the Miller Lab, Penn State.

Install from PyPi:
```
pip install tinscan
```

Clone and install from this repository:
```
git clone https://github.com/Adamtaranto/TE-insertion-scanner.git && cd TE-insertion-scanner && pip install -e .
```

### Example usage  

Find insertion events in genome A (target) relative to genome B (query). 


**Prepare Input Genomes**  


Split A and B genomes into two directories containing one scaffold per file.
Check that sequence names are unique within genomes.

```
tinscan-prep --adir data/A_target_split --bdir data/B_query_split\
-A data/A_target_genome.fasta -B data/B_query_genome.fasta 

```
Output: 
data/A_target_split/*.fa
data/B_query_split/*.fa


**Align Genomes**  


Align each scaffold from genome B onto each genome A scaffold.
Report alignments with >= 60% identity and length >= 100bp.

```
tinscan-align --adir data/A_target_split --bdir data/B_query_split \
--outdir A_Inserts --outfile A_Inserts_vs_B.tab \
--minIdt 60 --minLen 100 --hspthresh 3000

```
Output: 
A_Inserts/A_Inserts_vs_B.tab  


*Note:* Alignment tasks can be limited to a specified set of pairwise comparisons 
where appropriate (i.e. when homologous chromosome pairs are known between 
assemblies) using the option "--pairs". Comparisons are specified with a tab-delimited text file, where 
column 1 contains sequence names from genome A, and column 2 contains sequences 
from genome B.  

In the example *Chromosome_pairs.txt*, Chr A2 has been assembled 
as two scaffolds (B2,B3) in genome B.

```
#Chromosome_pairs.txt
A1	B1
A2	B2
A2	B3
A3	B4
```

**Find Insertions**  

Scan alignments for insertion events and report as GFF annotation of Genome A

```
tinscan-find --infile A_Inserts/A_Inserts_vs_B.tab \
--outdir A_Inserts --gffOut A_Inserts_vs_B_l100_id80.gff3 \
--maxInsert 50000 --minIdent 80 --maxIdentDiff 20

```
Output: 
A_Inserts/A_Inserts_vs_B_l100_id80.gff3

### Standard options

#### tinscan-prep

```
tinscan-prep --help
usage: tinscan-prep [-h] -A TARGET -B QUERY [--adir ADIR] [--bdir BDIR]
                    [-d OUTDIR]

Split multifasta genome files into directories for A and B genomes.

optional arguments:
  -h, --help            show this help message and exit
  -A TARGET, --target TARGET
                        Multifasta containing A genome.
  -B QUERY, --query QUERY
                        Multifasta containing B genome.
  --adir ADIR           A genome sub-directory within outdir
  --bdir BDIR           B genome sub-directory within outdir
  -d OUTDIR, --outdir OUTDIR
                        Write split directories within this directory.
                        (Default: cwd)
```

#### tinscan-align

```
tinscan-align --help
usage: tinscan-align [-h] --adir ADIR --bdir BDIR [--pairs PAIRS] [-d OUTDIR]
                     [--outfile OUTFILE] [--verbose] [--lzpath LZPATH]
                     [--minIdt MINIDT] [--minLen MINLEN]
                     [--hspthresh HSPTHRESH]

Align B genome (query) sequences onto A genome (target) using LASTZ.

optional arguments:
  -h, --help            show this help message and exit
  --adir ADIR           Name of directory containing sequences from A genome.
  --bdir BDIR           Name of directory containing sequences from B genome.
  --pairs PAIRS         Optional: Tab-delimited 2-col file specifying
                        target:query sequence pairs to be aligned
  -d OUTDIR, --outdir OUTDIR
                        Write output files to this directory. (Default: cwd)
  --outfile OUTFILE     Name of alignment result file.
  --verbose             If set report LASTZ progress.
  --lzpath LZPATH       Custom path to LASTZ executable if not in $PATH.
  --minIdt MINIDT       Minimum alignment identity to report.
  --minLen MINLEN       Minimum alignment length to report.
  --hspthresh HSPTHRESH
                        LASTZ min HSP threshold. Increase for stricter
                        matches.
```

#### tinscan-find

```
tinscan-find --help
usage: tinscan-find [-h] -i INFILE [--outdir OUTDIR] [--gffOut GFFOUT]
                    [--noflanks] [--maxTSD MAXTSD] [--maxInsert MAXINSERT]
                    [--minInsert MININSERT] [--qGap QGAP]
                    [--minIdent MINIDENT] [--maxIdentDiff MAXIDENTDIFF]

Parse whole genome alignments for signatures of transposon insertion.

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input file containing tab delimited LASTZ alignment
                        data.
  --outdir OUTDIR       Optional: Directory to write output to.
  --gffOut GFFOUT       Write features to this file as gff3.
  --noflanks            If set, do not report flanking hit regions in GFF.
  --maxTSD MAXTSD       Maximum overlap of insertion flanking sequences in
                        QUERY genome to be considered as target site
                        duplication. Flank pairs with greater overlaps will be
                        discarded Note: Setting this value too high may result
                        in tandem duplications in the target genome being
                        falsely classified as insertion events.
  --maxInsert MAXINSERT
                        Maximum length of sequence to consider as an insertion
                        event.
  --minInsert MININSERT
                        Minimum length of sequence to consider as an insertion
                        event. Note: If too short may detect small non-TE
                        indels.
  --qGap QGAP           Maximum gap allowed between aligned flanks in QUERY
                        sequence. Equivalent to target sequence deleted upon
                        insertion event.
  --minIdent MINIDENT   Minimum identity for a hit to be considered.
  --maxIdentDiff MAXIDENTDIFF
                        Maximum divergence in identity (to query) allowed
                        between insert flanking sequences.
```


# License

Software provided under MIT license.