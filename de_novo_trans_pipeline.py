### Adapted from https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
### Created by N.Vigneron, HES-SO Changins using Pycharm and patience ###
from functools import partial
import os
import glob
import subprocess
import multiprocessing
import time
import re
import sys
import pathlib as pl 
import argparse


###########################################################################################################################################

														##### Functions for shell #####
# Require in PATH variable trim_galore, bowtie2, bowtie2-build, run_rcorrector.pl, FilterUncorrectabledPEfastq.py

def run_rcorr(threads,forward, reverse):
	args = ["perl", "run_rcorrector.pl", "-t", threads ,"-1", forward, "-2",  reverse]
	r_corrected = subprocess.run(args)


def run_fix_read(threads, forward, reverse):
	filename = os.path.basename(forward)
	#mod_f = os.path.basename(forward)
	#pattern_SRA = re.compile("\w{3}\d{5,7}")
	#mod_f_2 = re.search(pattern_SRA, mod_f)
	#if mod_f_2: 
	#	filename = mod_f_2.group(0)
	args = ["python", "FilterUncorrectabledPEfastq.py","-1", forward, "-2",  reverse, "-s" , filename, "-t", threads]
	fix_read = subprocess.run(args)


def trim(threads, forward, reverse):
	args = ["trim_galore","--paired", "--phred33",  "--length", "36", "-q", "5", \
		"--stringency", "1", "-e", "0.1", "--path_to_cutadapt", "/usr/bin/cutadapt",  forward, reverse, "--gzip", "-j", "4"]
	trimed_read = subprocess.run(args)

def align(blacklist, threads, forward, reverse ):
	prefix = "_blacklist_"
	#mod_f= os.path.basename(forward)
	filename = os.path.basename(forward)
	#pattern_SRA = re.compile("\w{3}\d{5,7}")
	#mod_f_2 = re.search(pattern_SRA, mod_f)
	#if mod_f_2 :
	#	filename = mod_f_2.group(0)
	aligned_paired, aligned_unpaired = (filename + prefix + "Pv_Vv_paired_aligned" +".fq.gz"), (filename + prefix + "Pv_Vv_paired_unaligned"+".fq.gz")
	unaligned_paired, unaligned_unpaired = (filename + prefix + "Pv_Vv_unpaired_aligned"+".fq.gz"), (filename + prefix + "Pv_Vv_unpaired_unaligned" + ".fq.gz")
	args = ["bowtie2", "--quiet", "--very-sensitive-local", "--phred33 ", "-x", blacklist, "-1", forward, "-2", reverse, "--threads", threads,"--met-file", filename + "_bowtie2_metrics.txt", "--al-conc-gz", aligned_paired, "--un-conc-gz", unaligned_paired, "--al-gz", aligned_unpaired, "--un-gz", unaligned_unpaired ]
	align_read = subprocess.Popen(args)
	out = align_read.communicate()

def indexage(blacklist, threads):
	base = os.path.basename(blacklist)
	alias = os.path.splitext(base)[0]
	args = ["bowtie2-build" , "--large-index", blacklist, alias, "--threads", threads]
	hisat_align = subprocess.run(args)

														##### House-keeping functions #####

def cleaning(pattern):
	[os.remove(i) for i in glob.glob(pattern)]

def timer(start,end):
	hours, rem = divmod(end-start, 3600)
	minutes, seconds = divmod (rem, 60)
	print("Analysis took : {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

														#### Multiprocessing launcher functions ####

def multi_launch(function, threads, arg1, arg2, nb_tasks):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of cores to be used
		pool = multiprocessing.Pool(nb_tasks)
		intermediate = partial(function, str(threads))
		results = pool.starmap(intermediate, zip(arg1, arg2))
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

def multi_launch_bowtie(function, blacklist, threads, arg1, arg2, nb_tasks):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of cores to be used
		pool = multiprocessing.Pool(nb_tasks)
		intermediate = partial(function, blacklist, str(threads))
		results = pool.starmap(intermediate, zip(arg1, arg2))
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

###########################################################################################################################################
##### Beginning of instructions ######
usage = """
-dir working directory where directory and analysis will end up
-black path to file containing gene to exclude compose of SILVA rRNA, PN40024 12X v2 , P.viticola PV221
-t number of cpu use for the whole analysis, will be divided for parallelization speeding up analysis
-r directory where raw reads are stored default in exp_dir/sra if reads stores elsewhere use -r 
-job number of job running in parallel default = 10
"""

#### Get arguments from command line ####

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-dir", "--wrk_dir", dest="exp_dir", default="/path/to/De_novo", type =pl.Path, help="working directory where directory and analysis will end up")
parser.add_argument("-black", "--black_list", dest="black_list", default="/path/to/concatenate_SILVA_138.1_LSUParc_SSUParc_Pv_MBPM_Vv_12x_Cost.fasta", type =pl.Path, help="file containing black list")
parser.add_argument("-r", "--reads", dest="reads", default=None , type = pl.Path, help="directory where raw reads are stored")
parser.add_argument("-t", "--threads", dest="nb_threads", default= 6, type = int, help="Number of cpu use per job, Default = 6cpu")
parser.add_argument("-job", "--nb_job", dest="nb_job", default= 6, type = int, help="Number of job running in parallel, Default = 4 tasks in parallel")
args = parser.parse_args()

														#### Start of Analysis ####

print("Welcome to the De novo transcriptomic assembly pipeline\nPlease enter the parent directory you wish to use\n")

# Number of parallel tasks for the whole script #
number_of_tasks = args.nb_job
# Increase number of tasks for trim galore #
nb_threads = args.nb_threads

															# Set current directory #

#folder_path = input()
folder_path = str(args.exp_dir)
os.chdir(folder_path)

# Set an option to check is subdirectory of directory empty #
list_dir = os.listdir(os.getcwd())
to_remove = list_dir.remove("bowtied")
to_remove = list_dir.remove("sra")
for i in list_dir:
	if len(os.listdir(i)) == 0:
		print(f"Directory {i} empty")
	else:
		print(f"Caution directory {i} non empty, will interfere with analysis\n Please (re)move files in directory")
		ans = input("Do you want to continue ? press y:")
		if ans == "y":
			pass
		else:
			break

														### Checking blacklist ###

# Before running Hisat2 check if blacklist indexed #
#blacklist_path = input("Before launching analysis point to blacklist:\n")
blacklist_path = str(args.black_list)
if blacklist_path.endswith(".fasta"):
	pass
else:
	print("Path to blacklist incorrect\n")
	blacklist_path = input("Please enter correct path to blacklist")

# Get Database directory #
os.chdir(os.path.dirname(blacklist_path))

# Check if blacklist is indexed by checking number of files #
blacklist_name = os.path.splitext(blacklist_path)[0] + ".*"
test_blacklist = glob.glob(blacklist_name)
if len(test_blacklist) < 2:
	print("blacklist not indexed\n Indexing blacklist")
	indexage(blacklist_path, nb_threads*number_of_tasks)
else:
	print("blacklist already indexed, can move forward \n")


###########################################################################################################################################

															# Set directory #

# Create required directory (directory won't be erased if already exists) #
os.chdir(folder_path)
os.makedirs("rCorrected",exist_ok=True)
os.makedirs("fixed_read", exist_ok=True)
os.makedirs("trimmed_read", exist_ok=True)
os.makedirs("bowtied_read", exist_ok = True)

															# Launch rCorrector #

# Move to sra directory #
if args.reads is None:
	os.chdir("./sra")
else:
	os.chdir(str(args.reads))

# Make list of forward and reverse read files (required full path) #
list_1 = glob.glob(os.path.join(os.getcwd(),"*_R1.*"))
list_2 = glob.glob(os.path.join(os.getcwd(),"*_R2.*"))
list_1.sort(), list_2.sort()

# Change directory for further analysis #
os.chdir(os.path.join(folder_path, "rCorrected"))

# Multiprocessing launcher #
multi_launch(run_rcorr, nb_threads, list_1, list_2, number_of_tasks)

print("rCorrected analysis done !\n")

															# Launch Fixing read #

# Make list of forward and reverse read files (required full path) #
list_1 = glob.glob(os.path.join(os.getcwd(),"*_R1.*"))
list_2 = glob.glob(os.path.join(os.getcwd(),"*_R2.*"))
list_1.sort(), list_2.sort()
if len(list_1) == 0:
	sys.exit()

# Change directory for further analysis #
os.chdir("../fixed_read")

# Multiprocessing launcher #
multi_launch(run_fix_read,nb_threads, list_1, list_2, int(number_of_tasks/2))
print("Fixing read done !\n ")


															# Launch Trim Galore #

# Make list of forward and reverse read files (required full path) #
list_1 = glob.glob(os.path.join(os.getcwd(),"*_R1*"))
list_2 = glob.glob(os.path.join(os.getcwd(),"*_R2*"))
list_1.sort(), list_2.sort()
if len(list_1) == 0:
	sys.exit()
# Removing data from the rCorrected directory #
print("Cleaning")
os.chdir("../rCorrected")
cleaning("*.cor.*")

# Change directory for further analysis #
os.chdir("../trimmed_read")
# Multiprocessing launcher #
# Error if more than 4 threads used by tasks...
multi_launch(trim, int(nb_threads/2), list_1, list_2, number_of_tasks)


print("Trimming reads done !\n")

															# Launch Bowtie2 #

# Make list of forward and reverse read files (required full path) #
list_1 = glob.glob(os.path.join(os.getcwd(),"*1.cor_val_1.*"))
list_2 = glob.glob(os.path.join(os.getcwd(),"*2.cor_val_2.*"))
list_1.sort(), list_2.sort()
if len(list_1) == 0:
	sys.exit()
# Removing data from the fixed_read directory #
print("Cleaning")
os.chdir("../fixed_read")
cleaning("*.cor.fq.gz")

# Change directory for further analysis #
os.chdir("../bowtied_read")
# Multiprocessing launcher #
multi_launch_bowtie(align, os.path.splitext(blacklist_path)[0], nb_threads, list_1, list_2, number_of_tasks)

# Removing data from the trimmed_read directory #
print("Cleaning")
os.chdir("../trimmed_read")
cleaning("*.fq.gz")



print("End of script\nTransfer data to Trinity")
