#!/usr/bin/python

# PhyloCNV - estimating the abundance, gene-content, and phylogeny of microbes from metagenomes
# Copyright (C) 2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '0.0.2'

# Libraries
# ---------
import sys
import os
import subprocess
import resource
import gzip
from time import time
from platform import system

# Functions
# ---------

def add_paths(args):
	""" Add paths to external files and binaries """
	if system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % system())
	else:
		main_dir = os.path.dirname(os.path.abspath(__file__))
		args['bowtie2-build'] = '/'.join([main_dir, 'bin', system(), 'bowtie2-build'])
		args['bowtie2'] = '/'.join([main_dir, 'bin', system(), 'bowtie2'])
		args['samtools'] = '/'.join([main_dir, 'bin', system(), 'samtools'])
		args['pid_cutoffs'] = '/'.join([main_dir, 'data', 'pid_cutoffs.txt'])
		args['bad_gcs'] = '/'.join([main_dir, 'data', 'bad_cluster_ids.txt'])
		args['filter_bam'] = '/'.join([main_dir, 'filter_bam.py'])
		args['db'] = '/'.join([os.path.dirname(main_dir), 'ref_db/genome_clusters'])

def auto_detect_file_type(inpath):
	""" Detect file type [fasta or fastq] of <p_reads> """
	for line in iopen(inpath):
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % inpath)

def iopen(inpath):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzip.open(inpath)
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath))
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)

def genome_align(args):
	""" Use Bowtie2 to map reads to representative genomes from each genome cluster
	"""
	# Build command
	#	bowtie2
	command = '%s --no-unal ' % args['bowtie2']
	#   index
	command += '-x %s ' % '/'.join([args['out'], 'db', 'genomes'])
	#   specify reads
	if args['reads']: command += '-u %s ' % args['reads']
	#   speed/sensitivity
	command += '--%s ' % args['speed']
	#   threads
	command += '--threads %s ' % args['threads']
	#   file type
	if args['file_type'] == 'fasta': command += '-f '
	else: command += '-q '
	#   input file
	if (args['m1'] and args['m2']): command += '-1 %s -2 %s ' % (args['m1'], args['m2'])
	else: command += '-U %s' % args['m1']
	#   convert to bam
	command += '| %s view -b - ' % args['samtools']
	#   sort bam
	command += '| %s sort -f - %s ' % (args['samtools'], os.path.join(args['out'], 'genomes.bam'))
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def pileup(args):
	""" Filter alignments by % id, use samtools to create pileup, filter low quality bases, and write results to VCF file """
	# Build command
	#   percent id filtering
	command  = 'python %s %s %s %s | ' % (args['filter_bam'], '%s/genomes.bam' % args['out'], '/dev/stdout', args['mapid'])
	#   mpileup
	command += '%s mpileup -uv -A -d 10000 --skip-indels -B ' % args['samtools']
	#   quality filtering
	command += '-q %s -Q %s ' % (args['mapq'], args['baseq'])
	#   reference fna file
	command += '-f %s ' % ('%s/db/genomes.fa' % args['out'])
	#   input bam file
	command += '- '
	#   output vcf file
	command += '> %s ' % ('%s/genomes.vcf' % args['out'])
	# Run command
	if args['debug']: print("  running: %s") % command
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def read_ref_to_cluster(args, type):
	""" Read in map of scaffold id to genome-cluster id """
	ref_to_cluster = {}
	for line in open('/'.join([args['out'], 'db/%s.map' % type])):
		ref_id, cluster_id = line.rstrip().split()
		ref_to_cluster[ref_id] = cluster_id
	return ref_to_cluster

def split_vcf(args):
	""" Format vcf output for easy parsing """
	ref_to_cluster = read_ref_to_cluster(args, 'genomes')
	# open outfiles for each cluster_id
	outdir = '/'.join([args['out'], 'vcf'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	outfiles = {}
	for cluster_id in set(ref_to_cluster.values()):
		outfiles[cluster_id] = open('/'.join([outdir, '%s.vcf' % cluster_id]), 'w')
	# parse vcf into temorary vcf files for each cluster_id
	for line in open('/'.join([args['out'], 'genomes.vcf'])):
		if line[0] == '#': continue
		cluster_id = ref_to_cluster[line.split()[0]]
		outfiles[cluster_id].write(line)
	# close outfiles
	for file in outfiles.values():
		file.close()

def read_ref_bases(args, cluster_id):
	""" Read in reference genome by position """
	import Bio.SeqIO
	ref = []
	centroid_path = '/'.join([args['db'],cluster_id,'representative.fna.gz'])
	infile = gzip.open(centroid_path)
	for rec in Bio.SeqIO.parse(infile, 'fasta'):
		for pos in range(1, len(rec.seq)+1):
			ref.append([rec.id, pos, rec.seq[pos-1].upper()])
	return sorted(ref)

def write_snp_header(outfile):
	""" Write header for formatted SNP file """
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	outfile.write('\t'.join(fields)+'\n')

def write_snp_record(outfile, snp, ref):
	""" Write record for formatted SNP file """
	if ref:
		snp = {'ref_id': ref[0], 'ref_pos': str(ref[1]), 'ref_allele': ref[2],
		       'alt_allele': 'NA', 'cons_allele': 'NA', 'count_alleles': 1,
		       'depth': 0, 'count_ref': 0, 'count_alt': 0, 'ref_freq': 'NA'}
	fields = ['ref_id', 'ref_pos', 'ref_allele', 'alt_allele', 'cons_allele',
			  'count_alleles', 'count_ref', 'count_alt', 'depth', 'ref_freq']
	record = [str(snp[field]) for field in fields]
	outfile.write('\t'.join(record)+'\n')

def format_vcf(args):
	""" Format vcf files to snp files and fill in missing positions """
	outdir = '/'.join([args['out'], 'snps'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	for cluster_id in set(read_ref_to_cluster(args, 'genomes').values()):
		# open outfile
		outfile = gzip.open('/'.join([outdir, '%s.snps.gz' % cluster_id]), 'w')
		write_snp_header(outfile)
		# read sorted reference
		ref = read_ref_bases(args, cluster_id)
		ref_index = 0
		ref_length = len(ref)
		# write formatted records
		vcf_path = '/'.join([args['out'], 'vcf', '%s.vcf' % cluster_id])
		for snp in parse_vcf(vcf_path): # loop over formatted records from vcf
			snp_pos = [snp['ref_id'], int(snp['ref_pos'])]
			while snp_pos != ref[ref_index][0:2]: # fill in missing snp positions
				write_snp_record(outfile, None, ref[ref_index]) # write missing record
				ref_index += 1
			write_snp_record(outfile, snp, None) # write present record
			ref_index += 1
		while ref_index < ref_length: # fill in trailing snps
			write_snp_record(outfile, None, ref[ref_index]) # write trailing record
			ref_index += 1

def parse_vcf(inpath):
	""" Yields formatted records from VCF output """
	infile = open(inpath)
	for line in infile:
		r = line.rstrip().split()
		# get alt alleles
		alt_alleles = r[4].split(',')
		if '<X>' in alt_alleles: alt_alleles.remove('<X>')
		count_alleles = 1 + len(alt_alleles)
		# get allele counts
		info = dict([(_.split('=')) for _ in r[7].split(';')])
		counts = [int(_) for _ in info['I16'].split(',')[0:4]]
		# get consensus allele
		# *note: occassionally there are counts for alternate alleles, but no listed alternate alleles
		if sum(counts) == 0:
			cons_allele = 'NA'
		elif sum(counts[0:2]) >= sum(counts[2:4]):
			cons_allele = r[3]
		elif len(alt_alleles) == 0:
			cons_allele = 'NA'
		else:
			cons_allele = alt_alleles[0]
		# yield formatted record
		yield {'ref_id':r[0],
			   'ref_pos':r[1],
			   'ref_allele':r[3],
			   'count_alleles':count_alleles,
			   'alt_allele':alt_alleles[0] if count_alleles > 1 else 'NA',
			   'depth':sum(counts),
			   'count_ref':sum(counts[0:2]),
			   'count_alt':sum(counts[2:4]),
			   'cons_allele':cons_allele,
			   'ref_freq':'NA' if sum(counts) == 0 else sum(counts[0:2])/float(sum(counts))
			   }

def parse_file(inpath):
	""" Yields formatted records from SNPs output """
	infile = gzip.open(inpath)
	fields = next(infile).rstrip().split()
	for line in infile:
		yield dict([(i,j) for i,j in zip(fields, line.rstrip().split())])

def snps_summary(args):
	""" Get summary of mapping statistics """
	# store stats
	stats = {}
	for cluster_id in set(read_ref_to_cluster(args, 'genomes').values()):
		genome_length, covered_bases, total_depth, identity, maf = [0,0,0,0,0]
		for r in parse_file('/'.join([args['out'], 'snps/%s.snps.gz' % cluster_id])):
			genome_length += 1
			depth = int(r['depth'])
			if depth > 0:
				covered_bases += 1
				total_depth += depth
				if r['ref_allele'] == r['cons_allele']:
					identity += 1
				ref_freq = float(r['ref_freq'])
				maf += ref_freq if ref_freq <= 0.5 else 1 - ref_freq
		stats[cluster_id] = {'genome_length':genome_length,
							 'fraction_covered':covered_bases/float(genome_length),
							 'average_depth':total_depth/float(covered_bases) if covered_bases > 0 else 0,
							 'average_identity':identity/float(covered_bases) if covered_bases > 0 else 0,
							 'average_maf':maf/float(covered_bases) if covered_bases > 0 else 0
							 }
	# write stats
	fields = ['genome_length', 'fraction_covered', 'average_depth', 'average_identity', 'average_maf']
	outfile = open('/'.join([args['out'], 'snps_summary_stats.txt']), 'w')
	outfile.write('\t'.join(['cluster_id'] + fields)+'\n')
	for cluster_id in stats:
		record = [cluster_id] + [str(stats[cluster_id][field]) for field in fields]
		outfile.write('\t'.join(record)+'\n')

def max_mem_usage():
	""" Return max mem usage (Gb) of self and child processes """
	max_mem_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	max_mem_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
	return round((max_mem_self + max_mem_child)/float(1e6), 2)

def fetch_centroid(args, cluster_id):
	""" Get the genome_id corresponding to cluster centroid """
	inpath = '/'.join([args['db'], cluster_id, 'genomes.txt.gz'])
	infile = gzip.open(inpath)
	for line in infile:
		if line.split()[2] == 'Y':
			return line.split()[1]

def build_genome_db(args, genome_clusters):
	""" Build FASTA and BT2 database from genome cluster centroids """
	# fasta database
	outdir = '/'.join([args['out'], 'db'])
	if not os.path.isdir(outdir): os.mkdir(outdir)
	genomes_fasta = open('/'.join([args['out'], 'db', 'genomes.fa']), 'w')
	genomes_map = open('/'.join([args['out'], 'db', 'genomes.map']), 'w')
	db_stats = {'total_length':0, 'total_seqs':0, 'genome_clusters':0}
	for cluster_id in genome_clusters:
		if args['tax_mask'] and fetch_centroid(args, cluster_id) in args['tax_mask']:
			continue
		db_stats['genome_clusters'] += 1
		inpath = '/'.join([args['db'], cluster_id, 'representative.fna.gz'])
		infile = gzip.open(inpath)
		for line in infile:
			genomes_fasta.write(line)
			db_stats['total_length'] += len(line.rstrip())
			if line[0] == '>':
				sid = line.rstrip().lstrip('>').split()[0]
				genomes_map.write(sid+'\t'+cluster_id+'\n')
				db_stats['total_seqs'] += 1
	# print out database stats
	if args['verbose']:
		print("  total genomes: %s" % db_stats['genome_clusters'])
		print("  total contigs: %s" % db_stats['total_seqs'])
		print("  total base-pairs: %s" % db_stats['total_length'])
	# bowtie2 database
	inpath = '/'.join([args['out'], 'db', 'genomes.fa'])
	outpath = '/'.join([args['out'], 'db', 'genomes'])
	command = ' '.join([args['bowtie2-build'], inpath, outpath])
	if args['debug']: print("  running: %s" % command)
	process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = process.communicate()

def remove_tmp(args):
	""" Remove specified tmp files """
	if 'bowtie2_db' in args['remove']:
		os.remove('/'.join([args['out'], 'db/genomes.fa']))
		os.remove('/'.join([args['out'], 'db/genomes.fa.fai']))
		for file in os.listdir('%s/db' % args['out']):
			if file.split('.')[-1] == 'bt2' and file.split('.')[0] == 'genomes':
				os.remove('%s/db/%s' % (args['out'], file))
	if 'bam'in args['remove']:
		os.remove('%s/genomes.bam' % args['out'])
	if 'vcf'in args['remove']:
		import shutil
		shutil.rmtree('%s/vcf' % args['out'])
		os.remove('%s/genomes.vcf' % args['out'])

def run_pipeline(args):
	""" Run entire pipeline """
	
	add_paths(args) # Add paths to external files and binaries
	
	# Build genome database for selected GCs
	if args['build_db']:
		import species
		if args['verbose']: print("\nBuilding database of representative genomes")
		start = time()
		genome_clusters = species.select_genome_clusters(args)
		build_genome_db(args, genome_clusters)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# Use bowtie2 to map reads to a representative genome for each genome-cluster
	if args['align']:
		args['file_type'] = auto_detect_file_type(args['m1'])
		if args['verbose']: print("\nMapping reads to representative genomes")
		start = time()
		genome_align(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# Use mpileup to identify SNPs
	if args['pileup']:
		start = time()
		if args['verbose']: print("\nRunning mpileup")
		pileup(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# Split vcf into files for each GC, format, and report summary statistics
	if args['call']:
		start = time()
		if args['verbose']: print("\nFormatting output")
		split_vcf(args)
		format_vcf(args)
		snps_summary(args)
		if args['verbose']:
			print("  %s minutes" % round((time() - start)/60, 2) )
			print("  %s Gb maximum memory") % max_mem_usage()

	# Optionally remove temporary files
	if args['remove']:
		 remove_tmp(args)
