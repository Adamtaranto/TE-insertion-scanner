# Tinscan: TE-Insertion-Scanner

Scan whole genome alignments for transposon insertion signatures.

# Table of contents

* [Algorithm overview](#algorithm-overview)
* [Options and usage](#options-and-usage)
    * [Installing Tinscan](#installing-tinscan)
    * [Example usage](#example-usage)

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
  * Biopython

You can set up a conda environment with the required dependencies using the YAML files in this repo:

For ARM64 (Apple Silicon Macs) create a virtual intel env.

```bash
# For ARM64 Macs only
conda env create -f env_osx64.yml
conda activate tinscan-osx64
```


For all other operating systems use `environment.yml`

```bash
conda env create -f environment.yml
conda activate tinscan
```

With the conda env active you can now install tinscan.

1) Install from PyPi.
```bash
pip install tinscan
```

2) For the latest development version, clone and install from this repository.

```bash
git clone https://github.com/Adamtaranto/TE-insertion-scanner.git && cd TE-insertion-scanner && pip install -e ".[tests]"
```

### Example usage  

Find insertion events in genome A (target) relative to genome B (query). 


**Prepare Input Genomes**  


Split A and B genomes into two directories containing one scaffold per file.
Check that sequence names are unique within genomes.

```bash
tinscan-prep --adir data/A_target_split --bdir data/B_query_split\
-A data/A_target_genome.fasta -B data/B_query_genome.fasta 

```

Output: 
data/A_target_split/*.fa
data/B_query_split/*.fa


**Align Genomes**  


Align each scaffold from genome B onto each genome A scaffold.
Report alignments with >= 60% identity and length >= 100bp.

```bash
tinscan-align --adir data/A_target_split --bdir data/B_query_split \
--outdir A_Inserts --outfile A_Inserts_vs_B.tab \
--minIdt 60 --minLen 100 --hspthresh 3000

```

Output: 
A_Inserts/A_Inserts_vs_B.tab  


*Note:* Alignment tasks can be limited to a specified set of pairwise comparisons 
where appropriate (i.e. when homologous chromosome pairs are known between 
assemblies) using the option `--pairs`. 

Comparisons are specified with a 
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

```bash
tinscan-find --infile A_Inserts/A_Inserts_vs_B.tab \
--outdir A_Inserts --gffOut A_Inserts_vs_B_l100_id80.gff3 \
--maxInsert 50000 --minIdent 80 --maxIdentDiff 20

```

Output: 
A_Inserts/A_Inserts_vs_B_l100_id80.gff3

# License

Software provided under MIT license.