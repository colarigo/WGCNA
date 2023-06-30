### Created by N.Vigneron, HES-SO Changins using Pycharm and patience ###
import subprocess
import os, glob
import pathlib as pl
import sys
import argparse

# Require Trinity, jellyfish, salmon, bowtie2 in PATH variable

# Due to incompability (?) required to create bash script and execute it
def run_trinity(threads, forward, reverse):
	bash_script = "./launch_trinity.sh"
	#Creat bash script
	with open(bash_script, 'w') as handle:
		forward = ",".join(forward)
		reverse = ",".join(reverse)
		output = os.path.join(os.getcwd(),"trinity")
		handle.write("#!/bin/bash\n")
		handle.write(f"Trinity --no_version_check --seqType fq --no_normalize_reads --max_memory 50G --no_salmon --CPU {str(threads)} --left {forward} --right {reverse} --output {output}")
		#Make it executable
	subprocess.run(["chmod", "+x", bash_script])
	#Execute it
	subprocess.run([bash_script])

def run_supertranscripts(fasta):
	args = ["python", "/path/to/trinityrnaseq-v2.14.0.FULL/trinityrnaseq-v2.14.0/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py", "--incl_malign", "--trinity_fasta", fasta]
	run_sp = subprocess.run(args)

##### Beginning of instructions ######
usage = """
-dir working directory where directory and analysis will end up
-gen genome path
-gff path to annotation file in gff, gff3 or gtf format
-t number of cpu use for the whole analysis, will be divided for parallelization speeding up analysis
-r directory where raw reads are stored
-job number of job running in parallel default = 10
"""

#### Get arguments from command line ####

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-dir", "--wrk_dir", dest="exp_dir", default="/path/to/Output", type =pl.Path, help="working directory where directory and analysis will end up")
parser.add_argument("-r", "--reads", dest="reads", default=None , type = pl.Path, help="directory where align reads are stored", required = True)
parser.add_argument("-t", "--threads", dest="nb_threads", default= 20, type = int, help="Number of cpu use per job, Default = 20 cpu")
args = parser.parse_args()

os.chdir(str(args.exp_dir))
os.makedirs("trinity",exist_ok=True)
os.makedirs("./trinity/wrking_dir", exist_ok=True)
if len(args.reads) == None:
	os.chdir("./bowtied_read")
else:
	os.chdir(str(args.reads))
list_1 = glob.glob(os.path.join(os.getcwd(),"*_unpaired_aligned.fq.1.gz"))
list_2 = glob.glob(os.path.join(os.getcwd(),"*_unpaired_aligned.fq.2.gz"))
list_1.sort(), list_2 .sort()
if len(list_1) == 0 :
	sys.exit()
os.chdir("../trinity/wrking_dir")
run_trinity(args.nb_threads, list_1, list_2)

os.chdir("..")
list_1 = os.path.join(os.getcwd(),glob.glob("*Trinity.fasta$"))
if len(list_1) == 0 :
	sys.exit()
elif len(list_1) > 1:
	sys.exit()
os.chdir("..")
run_supertranscripts(list_1)