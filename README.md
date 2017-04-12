# Transposon insertion site prediction from genome alignments

Find alignment signatures characteristic of transposon insertion sites.

## Usage

Split target and query genomes into two directories containing one scaffold per file  
```bash
./SubsetFastaByName.py -i target.fa -d target_split --splitMode
./SubsetFastaByName.py -i query.fa -d query_split --splitMode
```  

Align query to target genome using LASTZ  
```bash
./LASTZ_genome_align.sh -z /usr/local/bin/lastz -t target_split -q query_split -i 60 -p query2target_90 -l 100 -d outdir
```  

Predict candidate insertion events and TSDs using insert finder
```bash
./insert-scanner.py -i outdir/query2target_90_concat.tab \ 
--maxTSD 100 \ 
--maxInsert 20000 \ 
--minInsert 100 \ 
--qGap 50 \ 
--minIdent 80 \ 
--maxIdentDiff 20 \ 
--outDir outdir \ 
--gffOut inserts.gff3
```