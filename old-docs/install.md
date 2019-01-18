### Requirements

* Operating systems: Unix, OSX
* Python >=2.7

Python modules (installed via setup.py):

* Numpy (>=1.7.0)
* BioPython (>=1.6.2)
* Pysam (>=0.8.1)
* Pandas (>=0.17.1)

External binaries (included with package):
* Bowtie2 (2.3.2)  
* samtools (1.4)  
* hs-blastn

### Speed and memory usage
Most steps of the program are quite fast, require < 2Gb of RAM, and are compatible with multithreading. For a typical gut metagenome, expect ~30 minutes for species profiling and 60-90 minutes each for calling nucleotide and gene copy number variants.

### Download the software

Clone the latest version from GitHub (recommended):  
`git clone https://github.com/snayfach/MIDAS`  
This makes it easy to update the software in the future using `git pull` as bugs are fixed and features are added

Or, download the latest [source code](https://github.com/snayfach/MIDAS/archive/master.zip), and unpack the project: `unzip MIDAS-master.zip`  

Or, download a stable release of the source code:  
[https://github.com/snayfach/MIDAS/releases](https://github.com/snayfach/MIDAS/releases)  

### Install dependencies

The dependencies can be installed by running the setup.py script:  
`python MIDAS/setup.py install` or  
`sudo python MIDAS/setup.py install` to install as a superuser, or  
`python MIDAS/setup.py install --user` to install locally  

Or, if you have pip installed, you can install these dependencies by running the following command:  
`pip install -U numpy biopython pysam pandas` or  
`sudo pip install -U numpy biopython pysam pandas` to install as a superuser, or  
`pip install --user numpy biopython pysam pandas` to install locally  

If you run into problems, check out these additional instructions:  
[Installing MIDAS on Mac OS X 10.11 and later](https://github.com/snayfach/MIDAS/issues/31)  
[Installing on CentOS/Ubuntu](install_other.md)

If you still have installation issues, please open a new [issue](https://github.com/snayfach/MIDAS/issues) on GitHub

### Update your environmental variables

This will enable you to run MIDAS scripts:  
`export PYTHONPATH=$PYTHONPATH:/path/to/MIDAS`  
`export PATH=$PATH:/path/to/MIDAS/scripts`   
`export MIDAS_DB=/path/to/midas_database` 

Be sure to replace '/path/to' with the appropriate path on your system  
Add these commands to your .bash_profile to avoid entering the commands in the future  

To learn how to download the default MIDAS database, click [here](ref_db.md)  
To learn how to build your own MIDAS database, click [here](build_db.md)   

### Testing

Run the following command to test MIDAS:  
`cd /path/to/MIDAS/test`  
`test_midas.py -vf`

This will take a few minutes to complete. If you don't have VSEARCH and HMMER installed, the last test (for building a custom database) will fail. This is OK and will not affect the rest of the program.  

If you do want to test the database building module, you'll need to add HMMER and VSEARCH installation directories to your PATH:  
`export PATH=$PATH:/vsearch/installation/directory`  
`export PATH=$PATH:/hmmer/installation/directory`

### Update MIDAS  
MIDAS changes over time as new features are added and bugs are fixed  

To update MIDAS, move to installation directory and pull the latest version:  
`cd MIDAS`  
`git pull` 

Or, download the latest stable release of the source code:  
[https://github.com/snayfach/MIDAS/releases](https://github.com/snayfach/MIDAS/releases)  
 
To be notified of these changes, click the "Watch" button in the upper right corner of the page, or sign up for the [mailing list](https://groups.google.com/forum/#!forum/midas-user-group) 

## Next step
[Step-by-step tutorial](docs/tutorial.md)
