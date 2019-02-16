try:
	from setuptools import setup
except:
	from distutils.core import setup
import os

setup(
	name = 'MIDAS',
	version = '0.0.0',
	description = 'An integrated pipeline and for estimating species and strain-level genomic variation from metagenomic data',
	license = 'GPL',
	author = 'Stephen Nayfach - witch changes by Boris Dimitrov',
	author_email='bdimitrov@chanzuckerberg.com',
	url='https://github.com/czbiohub/MIDAS-IGGdb',
	install_requires = ['biopython >= 1.62', 'numpy >= 1.7.0', 'pysam >= 0.8.1', 'pandas >= 0.17.1']
)
