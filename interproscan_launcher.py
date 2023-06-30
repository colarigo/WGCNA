### Created by N.Vigneron, HES-SO Changins using Pycharm and patience ###
import os, glob, os.path
import subprocess
import multiprocessing
import time
import re
import fileinput
import fnmatch
from Bio import SeqIO
import argparse

def chunk_itr(iterator, chunk_size):
    """
    Function to split fastq file into smallest files for faster processing
    From biopython solutions
    :param iterator: biopython FASTA file object
    :param chunk_size: number of sequences to write in each chunk
    """
    entry = True
    while entry:
        chunk = []
        while len(chunk) < chunk_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break
            chunk.append(entry)
        if chunk:
            yield chunk


def chunk_fasta(fasta_file, chunk_size, path_to_out):
    """
    Function for slicing a FASTA file in smaller FASTA files
    :param fasta_file: path to FASTA file to slice
    :param chunk_size: number of sequences in each FASTA chunks
    :param path_to_out: path to output directory
    :return: write chunks in output directory
    """
    rec_itr = SeqIO.parse(open(fasta_file), "fasta")
    file_name = os.path.splitext(os.path.basename(fasta_file))[0]
    for i, chunk in enumerate(chunk_itr(rec_itr, chunk_size)):
        out_chunk_name = os.path.join(path_to_out, "{0}_chunk{1}.fasta".format(file_name, i))
        with open(out_chunk_name, "w") as handle:
            SeqIO.write(chunk, handle, "fasta")

def interpros(fasta, threads):
    args = ["interproscan.sh", "-f", "HTML, TSV", "-cpu", threads, "-i", fasta, "-dp", "--goterms", "--pathways", "–iprlookup" ]
    interpro_annot = subprocess.Popen(args)
    out = interpro_annot.communicate()

def timer(start,end):
	hours, rem = divmod(end-start, 3600)
	minutes, seconds = divmod (rem, 60)
	print("Analysis took : {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

#Declaration of essential elements
chunk_size = 500
##### Beginning of instructions ######
usage = """
-i path to file containing protein sequence to annotate
-o directory where results are stored
-t number of cpu use for the whole analysis, will be divided for parallelization speeding up analysis
-job number of job running in parallel default = 4
"""
parser = argparse.ArgumentParser(usage=usage)

parser.add_argument("-i", "--input", dest="f_in", default="/path/to/fasta_file", type =pl.Path, help="path to file containing protein sequence to annotate")
parser.add_argument("-o", "--ouput", dest="f_out", default=None , type = pl.Path, help="directory where results are stored", required = True)
parser.add_argument("-t", "--threads", dest="nb_threads", default= 12, type = int, help="Number of cpu use per job, Default = 12 cpu")
parser.add_argument("-j", "--job", dest="nb_job", default= 4, type = int, help="Number job run in parallel, Default = 4")
args = parser.parse_args()

f_in = args.f_in
path_to_inter = args.f_out
nb_job = args.nb_job
nb_threads = args.nb_threads


#Get name of fasta file and create folder for interproscan input and output
f_in_name = os.path.basename(f_in).split(".")[0]
#Create parent directory with name of fasta files
os.makedirs(os.path.join(path_to_inter, f_in_name), exist_ok = True)
os.chdir(os.path.join(path_to_inter, f_in_name))
#Create subdirectories
os.makedirs(os.path.join(os.getcwd(), "interpro_in"), exist_ok = True)
pathout = os.path.join(os.getcwd(), "interpro_in")
os.makedirs(os.path.join(os.getcwd(), "interpro_out"), exist_ok = True)

start_1 = time.time()

# Remove * in fasta file
for line in fileinput.input(files = f_in, inplace=1, backup='.bak'):
	line = re.sub('\*$', '', line.rstrip())
	print(line)

chunk_fasta(f_in, chunk_size, pathout)


#Careful about filename extension,
#Might require modification
os.chdir(pathout)
chunk = glob.glob(os.path.join(os.getcwd(), "*.fasta"))
os.chdir("../interpro_out")

#Multiprocessing launcher
if __name__ == '__main__':
	#Specify number of cores to be used
	pool = multiprocessing.Pool(nb_job)
	results = pool.map(interpros, chunk, nb_threads)
	pool.close()
	pool.join()

end_1 = time.time()
timer(start_1,end_1)
print("Mission accomplished")

