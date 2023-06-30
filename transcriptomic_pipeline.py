### Created by N.Vigneron, HES-SO Changins using Pycharm and patience ###

from functools import partial
import os, glob, os.path
import subprocess
import multiprocessing
import time
import re
import pandas as pd
import pathlib as pl 
import argparse
import sys

### Require hisat2, samtools, featureCounts in PATH
### Functions declaration part ###
def indexage(genome,nb_threads):
	base = os.path.basename(genome)
	alias = os.path.splitext(base)[0]
	args = ["/path/to/hisat2-2.1.0/hisat2-build" , genome, alias, "-p", nb_threads]
	hisat_align = subprocess.run(args)

def run_trimo(nb_threads, forward, reverse):
	trim_path = "/usr/share/java/trimmomatic-0.36.jar"
	fname= os.path.basename(forward)
	#pattern_SRA = re.compile("\w{3}\d{7}")
	#mod_f_2 = re.search(pattern_SRA, mod_f)
	#if mod_f_2 : 
	#	fname = mod_f_2.group(0)
	args = ["java", "-jar", trim_path, "PE", forward, reverse, (fname + "_R1_paired.fq.gz"), (fname + "_R1_unpaired.fq.gz"), (fname + "_R2_paired.fq.gz"), (fname + "_R2_unpaired.fq.gz"), "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads",  "SLIDINGWINDOW:4:20", "LEADING:3", "TRAILING:3", "MINLEN:36", "-threads", nb_threads]
	r_trimo = subprocess.run(args)

def run_hisat2(genome,nb_threads, forward, reverse):
	fname= os.path.basename(forward)
	args = ["/path/to/hisat2-2.1.0/hisat2" ,"-x", genome, "-1", forward, "-2", reverse, "-S",  (fname + "_aligned.sam"), "--max-intronlen", "10000", "--dta", "-p", nb_threads]
	hisat_align = subprocess.run(args)

def sam_view(nb_threads, in_hisat):
	out_sam = os.path.splitext(os.path.basename(in_hisat))[0] + ".bam"
	args = ["samtools" ,"view", "-b", "-o", out_sam,  in_hisat, "--threads", nb_threads]
	sam_out = subprocess.run(args)
	print(f"Finished with {out_sam}")
	os.remove(in_hisat)


def sam_sort_name(nb_threads, in_bam):
	out_sam = os.path.splitext(os.path.basename(in_bam))[0] + "_sortbyName" +".bam"
	args = ["samtools" ,"sort", "-n", in_bam, "-o", out_sam, "--threads", nb_threads]
	sam_out = subprocess.run(args)
	print(f"Finished with {out_sam}")


def sam_fixmate(nb_threads, in_bam):
	out_sam = os.path.splitext(os.path.basename(in_bam))[0] + "_fixmate" +".bam"
	args = ["samtools" ,"fixmate", "-O", "BAM", "-m", in_bam, out_sam, "--threads", nb_threads]
	sam_out = subprocess.run(args)
	print(f"Finished with {out_sam}")

def sam_sort(nb_threads, in_bam):
	out_sam = os.path.splitext(os.path.basename(in_bam))[0] + "_sort" +".bam"
	args = ["samtools" ,"sort", in_bam, "-o", out_sam, "--threads", nb_threads]
	sam_out = subprocess.run(args)
	print(f"Finished with {out_sam}")

def sam_markdup(nb_threads, in_bam):
	out_sam = os.path.splitext(os.path.basename(in_bam))[0] + "_markdup" +".bam"
	args = ["samtools" ,"markdup", "-r", "-s", "-O", "BAM",  in_bam, out_sam, "--threads", nb_threads]
	sam_out = subprocess.run(args)
	print(f"Finished with {out_sam}")

def feat_count(threads, gff, sample):
	out_count = os.path.splitext(os.path.basename(sample))[0] + ".txt"
	args = ["featureCounts" ,"-a", gff, "-o", out_count, "-t", "gene", "-g", "ID", "-T", str(threads), "-p", "-B", sample]
	count_out = subprocess.run(args)

def timer(start,end):
	hours, rem = divmod(end-start, 3600)
	minutes, seconds = divmod (rem, 60)
	print("Analysis took : {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

def multi_launch_sam(function, arg1, nb_tasks, nb_threads):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of cores to be used
		pool = multiprocessing.Pool(nb_tasks)
		intermediate = partial(function, nb_threads) 
		results = pool.map(intermediate, arg1)
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

def multi_launch_trimo(function, arg1, arg2, nb_tasks, nb_threads):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of cores to be used
		pool = multiprocessing.Pool(nb_tasks)
		intermediate = partial(function, nb_threads)
		results = pool.starmap(intermediate, zip(arg1, arg2))
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

def multi_launch_hisat(function, genome, arg1, arg2, nb_tasks,  nb_threads):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of cores to be used
		pool = multiprocessing.Pool(nb_tasks)
		intermediate = partial(function, genome,  nb_threads)
		results = pool.starmap(intermediate, zip(arg1, arg2))
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

def cleaning(pattern):
	to_remove = glob.glob(pattern)
	for i in to_remove:
		os.remove(i)

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
parser.add_argument("-dir", "--wrk_dir", dest="exp_dir", default="/path/to/transcriptomic", type =pl.Path, help="working directory where directory and analysis will end up")
parser.add_argument("-gff", "--gen_annot", dest="gen_annot", default=None, type =pl.Path, help="gff file containing gene annotation", required = True)
parser.add_argument("-gen", "--gen_path", dest="gen_path", default=None, type = pl.Path, help="path to genome in fasta", required = True)
parser.add_argument("-r", "--reads", dest="reads", default=None , type = pl.Path, help="directory where raw reads are stored", required = True)
parser.add_argument("-t", "--threads", dest="nb_threads", default= 12, type = int, help="Number of cpu use per job, Default = 12 cpu")
parser.add_argument("-job", "--nb_job", dest="nb_job", default= 4, type = int, help="Number of job running in parallel, Default = 4 tasks in parallel")
args = parser.parse_args()

#### Beginning of analysis ####
folder_path = str(args.exp_dir)
path_to_gff= str(args.gen_annot)
genome_path = str(args.gen_path)
number_of_threads = args.nb_threads
number_of_tasks = args.nb_job
sra_path = str(args.reads)

print("Welcome to the transcriptomic pipeline\n")
#### Setting genome and gff path and verification on the presence of index and format of files ####
#Set working directory
#If folder displays is the correct one moving on otherwise it is possible to change it using input()
print(f"Current output directory is {folder_path}\n")
#Using while in case several mistakes have been made to enter a correct path
new_folder_path = input("If you wish to change press y otherwise press any key:")
if new_folder_path == "y":
	folder_path = input("Enter the new output directory:\n")
	while not folder_path.startswith("/") or folder_path.startswith('"'):
		folder_path = input("Path incorrect\nEnter the new output directory:\n")
else:
	pass

### Checking genome ###
#Before running Hisat2 check if genome indexed
if genome_path.endswith(".fasta"):
	pass
else:
	print("Path to genome incorrect\n")
	genome_path = input("Please enter correct path to genome")

#Get Database directory
os.chdir(os.path.dirname(genome_path))
#Check if genome is indexed by checking number of files
genome_name = os.path.splitext(genome_path)[0] + ".*"
test_genome = glob.glob(genome_name)
if len(test_genome) < 3:
	print("Genome not indexed\n Indexing genome")
	indexage(genome_path, str(number_of_threads*number_of_tasks))
else:
	print("Genome already indexed, can move forward \n")

### Checking gff file ###
#Get gff file
# Check if gff file exists and if it ends with gff3/gff/gtf
if os.path.isfile(path_to_gff) is False or not re.search("\.gff3|\.gff|\.gtf", path_to_gff):
	print("Gff file not found or the path does not direct toward an annotation file\nPlease enter the full and corrected path")
	path_to_gff = input()
else:
	print("Gff file exists, moving on")	
	pass

### Verification of genome and gff file done ###
start_1 = time.time()
#### Starting first part of analysis #####
###  Trimmomatic and Hisat2 ###
os.chdir(folder_path)
#Create required directory (directory won't be erased if already exists)
os.makedirs("trim_out",exist_ok=True)
os.makedirs("hisat2_out", exist_ok=True)
os.makedirs("sam_out",exist_ok = True)
os.makedirs("featurecounts", exist_ok = True)


#Make list of forward and reverse read files (required full path)
os.chdir(sra_path)
list_1 = glob.glob(os.path.join(sra_path, "*R1.fastq.gz"))
list_2 = glob.glob(os.path.join(sra_path, "*R2.fastq.gz"))
list_1.sort(), list_2.sort()
if len(list_1) == 0 :
	sys.exit()
##After getting path of sra files changes to trimmomatic output
os.chdir(os.path.join(folder_path,"trim_out"))
multi_launch_trimo(run_trimo, list_1, list_2, number_of_tasks, str(number_of_threads))

#Let's clean a bit
print("Trimmming done cleaning and moving on")
cleaning("*_unpaired.fq.gz")

#Get list of trim read files and change directory
list_1 = glob.glob(os.path.join(os.getcwd(), "*R1_paired.fq.gz"))
list_2 = glob.glob(os.path.join(os.getcwd(), "*R2_paired.fq.gz"))
list_1.sort(), list_2.sort()
os.chdir("../hisat2_out")
if len(list_1) == 0 :
	sys.exit()

#Start Hisat2 analysis
multi_launch_hisat(run_hisat2, os.path.splitext(genome_path)[0], list_1, list_2, number_of_tasks, str(number_of_threads))
os.chdir("../trim_out")
cleaning("*_paired.fq.gz")
os.chdir("../hisat2_out")
#
###### Using samtools to modify files #####
####   Started first sorting ###
##     Get all aligments path #
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned.sam"))
if len(list_1) == 0 :
	sys.exit()
os.chdir("../sam_out")
##Convert sam to bam
multi_launch_sam(sam_view, list_1, number_of_tasks, str(number_of_threads))
print("Samtools view analysis done\nStarting Samtools sorting")
##Get bam file
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned.bam"))
if len(list_1) == 0 :
	sys.exit()
##Samtools sort
multi_launch_sam(sam_sort_name, list_1, number_of_tasks, str(number_of_threads))
print("Samtools sort analysis done\nStarting Samtools fixmate and a bit of cleaning")
## Clean a bit
cleaning("*_aligned.bam")
os.chdir("../hisat2_out")
# cleaning("*_aligned.sam")
os.chdir("../sam_out")
#
#### Launch fixmate ###
##   Get sorted files #
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned_sortbyName.bam"))
if len(list_1) == 0 :
	sys.exit()
##Samtools fixamte
multi_launch_sam(sam_fixmate, list_1, number_of_tasks, str(number_of_threads))
print("Samtools fixmate analysis done\nStarting Samtools sort")
##Bit of cleaning again
cleaning("*_aligned_sortbyName.bam")
#
#### Second sorting ###
##   Get fixamte files #
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned_sortbyName_fixmate.bam"))
if len(list_1) == 0 :
	sys.exit()
##Samtools sort
multi_launch_sam(sam_sort, list_1, number_of_tasks, str(number_of_threads))
print("Files sorted by order\nStarting remove duplicate and a bit of cleaning")
##Bit of cleaning...
cleaning("*_aligned_sortbyName_fixmate.bam")
#
#### Launch markdup ###
##   Get fixmate sorted files #
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned_sortbyName_fixmate_sort.bam"))
if len(list_1) == 0 :
	sys.exit()
##Samtools markdup
multi_launch_sam(sam_markdup, list_1, number_of_tasks, str(number_of_threads))
##Bit of cleaning again and again
cleaning("*_aligned_sortbyName_fixmate_sort.bam")
print("Duplicates remove and cleaning done\nFinal sorting")
#
#### Last sorting ###
list_1 = glob.glob(os.path.join(os.getcwd(), "*_aligned_sortbyName_fixmate_sort_markdup.bam"))
if len(list_1) == 0 :
	sys.exit()
##Samtools sort by name
multi_launch_sam(sam_sort_name, list_1, number_of_tasks, str(number_of_threads))
##Last cleaning
cleaning("*_aligned_sortbyName_fixmate_sort_markdup.bam")
print("Garbage has been collected")
print("Files have been sorted, fixmate and duplicates were removed\nReady for further analysis")

print("Starting featureCounts")
samples = glob.glob(os.path.join(os.getcwd(), "*.bam"))
if len(list_1) == 0 :
	sys.exit()
os.chdir("../featurecounts")
for i in samples:
	feat_count(number_of_threads, path_to_gff, i)
for i in glob.glob("*.txt"):
	data = pd.read_csv(i, sep='\t', skiprows= [0],header= 0,usecols = [0,6], names= ['Gene','Counts'])
	names_df = os.path.splitext(i)[0] +'_counts.txt'
	data.to_csv(names_df, sep= '\t', index=False, header=False)
os.chdir("../sam_out")
cleaning("*.bam")

end_1 = time.time()
print("Pipeline finished\nEnjoy :)")
timer(start_1, end_1)



