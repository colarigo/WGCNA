import os
import subprocess
import multiprocessing
import time
from functools import partial
import pathlib as pl 
import argparse
import pandas as pd


def timer(start,end):
	hours, rem = divmod(end-start, 3600)
	minutes, seconds = divmod (rem, 60)
	print("Analysis took : {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

def clustering(opt,path):
	os.chdir(path)
	output = "Results" + "_" + opt + "_" + os.path.basename(path)
	if int(opt) != 1:
		args = ["python","path/to/clust-1.17.0/clust.py", path, "-n 4", "-t", opt, "-o", output]
	else:
		args = ["python","path/to/clust-1.17.0/clust.py", path, "-n 4", "-o", output]
	cluster = subprocess.run(args)
	print(f"Analysis done for {os.path.basename(path)}")

def multi_launch_partial(function, path, nb_tasks, opt):
	if __name__ == '__main__':
		#start timer
		start_1 = time.time()
		#Specify number of job to execute in parallel
		pool = multiprocessing.Pool(nb_tasks)
		int_genome = partial(function, str(opt))
		results = pool.map(int_genome, path)
		pool.close()
		pool.join()
		end_1 = time.time()
		timer(start_1,end_1)

##### Beginning of instructions ######
usage = """
-dir working directory
-opt tightness of cluster
-annot gene annotation
"""

#### Get arguments from command line ####
parser = argparse.ArgumentParser(usage=usage)
parser.add_argument("-dir", "--wrk_dir", dest="exp_dir", type =pl.Path, help="working directory where directory and analysis will end up",required = True)
parser.add_argument("-opt",dest="opt", default = 1, type =int, help="degree of module tightness, integer between 0 and 10")
parser.add_argument("-annot",dest="annot", type =pl.Path, help="path to gene annotation", required = True)
args = parser.parse_args()

#### Beginning of analysis ####
source = str(args.exp_dir)
tight = str(args.opt)
folder_to = []
for root, dirs, files in os.walk(source, topdown=False):
	for name in dirs:	
		if any(name.startswith(x) for x in ["Data", "Processed", "Input", "Results"]):
			pass
		else:
			folder_to.append(os.path.join(root, name))

multi_launch_partial(clustering, folder_to, len(folder_to), tight)


os.chdir(source)
cluster = []
for root, dirs, files in os.walk(os.getcwd(), topdown=False):
	for name in files:
		if name == "Clusters_Objects.tsv":
			# Remove other analysis with different opt int
			if os.path.basename(os.path.dirname(os.path.join(root,name))).startswith("".join(["Results_", str(tight)])) :
				cluster.append(os.path.join(root,name))


# list of dataframe
list_df = []
# Gene annotation
annot = pd.read_csv(args.annot, sep="\t")
for i in range(len(cluster)):
    if os.path.getsize(cluster[i]) > 0:
        df1 = pd.read_csv(cluster[i], sep="\t")
        # Remove first row
        df1 = df1.iloc[1: , :]
        # Rename columns [ x for x in y] equivalent of paste
        dirs = os.path.basename(os.path.dirname(os.path.dirname(cluster[i])))
        df1.columns = [dirs + "_" + str(x) for x in list(range(1, len(df1.columns)+ 1))]
        # Long format
        df1 = pd.melt(df1, var_name='clust', value_name='gene')
        # Remove Nan
        df1 = df1.dropna(subset=['gene'])
        # Add gene annotation and process
        inner_join = pd.merge(df1, annot, on ='gene', how ='inner')
        list_df.append(inner_join)
pd.concat(list_df).to_csv("".join(["Clusters_",os.path.basename(os.getcwd()),".txt"]), sep="\t", index= False)


