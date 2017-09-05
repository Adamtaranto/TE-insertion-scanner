# Tinscan: TE-Insertion-Scanner

Scan whole genome alignments for transposon insertion signatures.


# Table of contents

* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing Tinscan](#installing-tinscan)
    * [Example usage](#example-usage)
    * [Standard options](#standard-options)

# Algorithm overview

  1. 
  2.  


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
tinscan-prep --outdir GenA_Insertions -A genomeA.fa -B genomeB.fa
```
Output: 
GenA_Insertions/A_target_split/*.fa
GenA_Insertions/B_query_split/*.fa


**Align Genomes**  

Align each scaffold from genome B onto each genome A scaffold.
Report alignments with >= 60% identity and length >= 100bp.

```
tinscan-align -d GenA_Insertions --outfile alignB2A.tab \
-A GenA_Insertions/A_target_split -B GenA_Insertions/B_query_split \
-i 60 -p alignB2A -l 100

```
Output: 
/alignB2A.tab

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
tinscan-find -i GenA_Insertions/alignB2A.tab
--maxTSD 100 \
--maxInsert 20000 \
--minInsert 100 \
--qGap 50 \
--minIdent 80 \
--maxIdentDiff 20 \
--outDir GenA_Insertions \
--gffOut GenA_Insertions_vs_GenB.gff3

```

### Standard options




# License

Software provided under MIT license.