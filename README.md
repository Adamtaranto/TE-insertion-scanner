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
  4. At least one flanking alignment in A satisfies the threshold --minIdent, and differs from its mate by no more than --maxIdentDiff %
  5. Log flanks and candidate insertion.
  6. Attempt to infer TSDs from the internal overlap of flanking alignments in B genome.

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
assemblies) using the option "--pairs". Comparisons are specified with a 
tab-delimited text file, where column 1 contains sequence names from genome A, and 
column 2 contains sequences from genome B.  

In the example *Chromosome_pairs.txt*, Chr A2 has been assembled 
as two scaffolds (B2, B3) in genome B.

```
#Chromosome_pairs.txt
A1    B1
A2    B2
A2    B3
A3    B4
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
Usage: tinscan-prep [-h] -A TARGET -B QUERY [--adir ADIR] [--bdir BDIR]
                    [-d OUTDIR]

Split multifasta genome files into directories for A and B genomes.

Optional arguments:
  -h, --help        Show this help message and exit.
  -A , --target     Multifasta containing A genome.
  -B , --query      Multifasta containing B genome.
  --adir            A genome sub-directory within outdir
  --bdir            B genome sub-directory within outdir
  -d , --outdir     Write split directories within this directory.
                    (Default: cwd)
```

#### tinscan-align

```
tinscan-align --help
Usage: tinscan-align [-h] --adir ADIR --bdir BDIR [--pairs PAIRS] [-d OUTDIR]
                     [--outfile OUTFILE] [--verbose] [--lzpath LZPATH]
                     [--minIdt MINIDT] [--minLen MINLEN]
                     [--hspthresh HSPTHRESH]

Align B genome (query) sequences onto A genome (target) using LASTZ.

Optional arguments:
  -h, --help        Show this help message and exit
  --adir            Name of the directory containing sequences from A genome.
  --bdir            Name of the directory containing sequences from B genome.
  --pairs           Optional: Tab-delimited 2-col file specifying
                    target:query sequence pairs to be aligned
  -d , --outdir     Write output files to this directory. (Default: cwd)
  --outfile         Name of alignment result file.
  --verbose         If set report LASTZ progress.
  --lzpath          Custom path to LASTZ executable if not in $PATH.
  --minIdt          Minimum alignment identity to report.
  --minLen          Minimum alignment length to report.
  --hspthresh       LASTZ min HSP threshold. Increase for stricter matches.
```

#### tinscan-find

```
tinscan-find --help
Usage: tinscan-find [-h] -i INFILE [--outdir OUTDIR] [--gffOut GFFOUT]
                    [--noflanks] [--maxTSD MAXTSD] [--maxInsert MAXINSERT]
                    [--minInsert MININSERT] [--qGap QGAP]
                    [--minIdent MINIDENT] [--maxIdentDiff MAXIDENTDIFF]

Parse whole genome alignments for signatures of transposon insertion.

Optional arguments:
  -h, --help        Show this help message and exit.
  -i , --infile     An input file containing tab-delimited LASTZ alignment
                    data.
  --outdir          Optional: Directory to write output to.
  --gffOut          Write features to this file as GFF3.
  --noflanks        If set, do not report flanking hit regions in GFF.
  --maxTSD          Maximum overlap of insertion flanking sequences in
                    QUERY genome to be considered as target site
                    duplication. Flank pairs with greater overlaps will be
                    discarded Note: Setting this value too high may result
                    in tandem duplications in the target genome being
                    falsely classified as insertion events.
  --maxInsert       The maximum sequence length to consider as an insertion 
                    event.
  --minInsert       Minimum sequence length to consider as an insertion
                    event. Note: If too short may detect small non-TE
                    indels.
  --qGap            Maximum gap allowed between aligned flanks in QUERY
                    sequence. Equivalent to target sequence deleted upon
                    insertion event.
  --minIdent        Minimum identity for a hit to be considered.
  --maxIdentDiff    Maximum divergence in identity (to query) allowed
                    between insert flanking sequences.
```

# License

Software provided under MIT license.