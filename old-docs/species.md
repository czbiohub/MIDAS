## Metagenomic species profiling

This script will map metagenomic reads to a database of phylogenetic marker genes using HS-BLASTN. Mapped reads are used estimate the read depth and relative abundance of 5,952 bacterial species. Reads are mapped according to gene-specific, species-level mapping thresholds (94.5-98% DNA identity). Reads that map equally well to 2 or more species are probabalistically assigned.

This module is required for downstream modules (i.e. genes and snps) if you plan to automatically identify abundant species with no prior knowlege of the community. If you already have a list of species ids that you're interested, this step can be skipped.

## Usage
```
Usage: run_midas.py species outdir [options]

positional arguments:
  outdir             Path to directory to store results. 
                     Name should correspond to unique sample identifier.

optional arguments:
  -h, --help         show this help message and exit
  -1 M1              FASTA/FASTQ file containing 1st mate if using paired-end reads.
                     Otherwise FASTA/FASTQ containing unpaired reads.
                     Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
  -2 M2              FASTA/FASTQ file containing 2nd mate if using paired-end reads.
                     Can be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)
  -n MAX_READS       Number of reads to use from input file(s) (use all)
  -t THREADS         Number of threads to use for database search (1)
  -d DB              Path to reference database
                     By default, the MIDAS_DB environmental variable is used
  --remove_temp      Remove temporary files, including BLAST output.
                     Useful for reducing disk space of MIDAS output
  --word_size INT    Word size for BLAST search (28)
                     Use word sizes > 16 for greatest efficiency.
  --mapid FLOAT      Discard reads with alignment identity < MAPID
                     By default gene-specific species-level cutoffs are used
                     Values between 0-100 accepted
  --aln_cov FLOAT    Discard reads with alignment coverage < ALN_COV (0.75)
                     Values between 0-1 accepted
  --read_length INT  Trim reads to READ_LENGTH and discard reads with length < READ_LENGTH
                     By default, reads are not trimmed or filtered
```

## Examples
1) run with defaults using a paired-end metagenome:  
`run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz -2 /path/to/reads_2.fq.gz`

2) run using a single-end metagenome with 4 CPUs and only 4M reads:  
`run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz -t 4 -n 4000000`

3) run with exactly 80 base-pair reads:  
`run_midas.py species /path/to/outdir -1 /path/to/reads_1.fq.gz --read_length 80`

## Output files
<b>species_profile.txt:</b> tab-delimited with header; each line contains the abundance values for 1 species (5,952 total species) sorted by decreasing relative abundance. 

<b>log.txt:</b> log file containing parameters used

<b>temp:</b> directory of intermediate files; run with `--remove_temp` to remove these files
  
## Output formats
<b>species_profile.txt:</b>

  * species_id: species identifier
  * count_reads: number of reads mapped to marker genes
  * coverage: estimated genome-coverage (i.e. read-depth) of species in metagenome
  * relative_abundance: estimated relative abundance of species in metagenome
  
## Memory usage
* < 1.5 Gb for most samples

## Speed
* ~5,000 reads/second for 100-bp reads when using default parameters
* Use `-n` and `-t` to increase throughput
* Note than using `-n` will result in underestimates of species genome-coverage in the full metagenome
* We found that about 1 million reads was sufficient to precisely estimate species relative abundance for a gut community

## Next step
[Merge species abundance across samples](merge_species.md)
