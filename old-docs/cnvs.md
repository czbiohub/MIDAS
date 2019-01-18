## CAUTION:  This verbatim copy of the original [MIDAS tool](https://github.com/snayfach/MIDAS) documentation is OUT OF DATE for the MIDAS-IGGdb update, and only included as a temporary reference during construction of MIDAS-IGGdb.


## Metagenomic pan-genome profiling
This script will map metagenomic reads to bacterial pangenomes and quantify these genes in your data. You can either target one or more specific species (`--species_id`), or provide this script a species abundance file. 

The pipeline can be broken down into the following steps:
  
  * build a database of pangenomes for abundance bacterial species
  * map high-quality metagenomic reads to the database
  * use mapped reads to quantify pangenome genes

## Usage
```
Usage: run_midas.py genes <outdir> [options]

positional arguments:
  outdir                Path to directory to store results.
                        Directory name should correspond to sample identifier

optional arguments:
  -h, --help            show this help message and exit
  --remove_temp         Remove intermediate files generated by MIDAS (False).
                        Useful to reduce disk space of MIDAS output

Pipeline options (choose one or more; default=all):
  --build_db            Build bowtie2 database of pangenomes
  --align               Align reads to pangenome database
  --call_genes          Compute coverage of genes in pangenome database

Database options (if using --build_db):
  -d DB                 Path to reference database
                        By default, the MIDAS_DB environmental variable is used
  --species_cov FLOAT   Include species with >X coverage (3.0)
  --species_topn INT    Include top N most abundant species
  --species_id CHAR     Include specified species. Separate ids with a comma

Read alignment options (if using --align):
  -1 M1                 FASTA/FASTQ file containing 1st mate if using paired-end reads.
                        Otherwise FASTA/FASTQ containing unpaired reads.
                        Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
  -2 M2                 FASTA/FASTQ file containing 2nd mate if using paired-end reads.
                        Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
  --interleaved         FASTA/FASTQ file in -1 are paired and contain forward AND reverse reads
  -s {very-fast,fast,sensitive,very-sensitive}
                        Alignment speed/sensitivity (very-sensitive)
  -m {local,global}     Global/local read alignment (local)
  -n MAX_READS          # reads to use from input file(s) (use all)
  -t THREADS            Number of threads to use (1)

Quantify genes options (if using --call_genes):
  --readq INT           Discard reads with mean quality < READQ (20)
  --mapid FLOAT         Discard reads with alignment identity < MAPID (94.0)
  --aln_cov FLOAT       Discard reads with alignment coverage < ALN_COV (0.75)
  --trim INT            Trim N base-pairs from 3'/right end of read (0)
```

## Examples

1) run entire pipeline using defaults:  
`run_midas.py genes /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz`

2) run entire pipeline for a specific species:  
`run_midas.py genes /path/to/outdir --species_id Bacteroides_vulgatus_57955 -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz`

3) just align reads, use faster alignment, only use the first 10M reads, use 4 CPUs:  
`run_midas.py genes /path/to/outdir --align -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz -s very-fast -n 10000000 -t 4`

4) just quantify genes, keep reads with >=95% alignment identity and reads with an average quality-score >=30:  
`run_midas.py genes /path/to/outdir --call_genes --mapid 95 --readq 20`

## Output files
 
<b>output:</b> directory of per-species output files; files are tab-delimited, gzip-compressed, with header.  
<b>species.txt:</b> list of species_ids included in local database  
<b>summary.txt:</b> tab-delimited with header; summarizes alignment results per-species  
<b>log.txt:</b> log file containing parameters used  
<b>temp:</b> directory of intermediate files; run with `--remove_temp` to remove these files  
  
## Output formats  

<b>output/\<species\_id>.genes.gz</b>
  
  * gene_id: id of non-redundant gene used for read mapping; 'peg' and 'rna' indicate coding & RNA genes respectively  
  * count_reads: number of aligned reads to gene\_id after quality filtering  
  * coverage: average read-depth of gene_id based on aligned reads (# aligned bp / gene length in bp)  
  * copy_number: estimated copy-number of gene\_id based on aligned reads (coverage of gene\_id / median coverage of 15 universal single copy genes)  

<b>summary.txt</b>
  
  * species_id: species id  
  * pangenome_size: number of non-redundant genes in reference pan-genome  
  * covered_genes: number of genes with at least 1 mapped read  
  * fraction_covered: proportion of genes with at least 1 mapped read  
  * mean_coverage: average read-depth across genes with at least 1 mapped read  
  * marker_coverage: median read-depth across 15 universal single copy genes  
  * aligned_reads: number of aligned reads BEFORE quality filtering  
  * mapped_reads: number of aligned reads AFTER quality filtering  
  
## Memory usage  
* Memory usage will depend on the number of species you search and the number of reference genomes sequenced per species.
* In practice, peak memory usage will not exceed 1 Gb for most samples

## Speed
* Speed will depend on the number of species you search and the number of reference genomes sequenced per species. 
* For a single species with 1 reference genome, expect ~16,000 reads/second
* Use `-n` and `-t` to increase throughput

## Next step
[Merge results across samples](merge_cnvs.md)
