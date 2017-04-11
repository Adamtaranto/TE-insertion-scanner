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
./LASTZ_genome_align.sh target_split query_split 90 query2target_90
```  

Predict candidate insertion events and TSDs using insert finder
```bash
./insert-finder.py -i query2target_90.tab
```