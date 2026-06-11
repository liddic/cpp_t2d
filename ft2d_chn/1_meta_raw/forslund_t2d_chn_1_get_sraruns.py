# # # # # # # # # # # # # # # #
# # Download SRA fastq files - Forslund T2D Chinese cohort  metagenomics
# # # # # # # # # # # # # # # #
#$ sinfo -N -o "%N %C %e %m"   # get cluster usage info
#$ cd <check run folder!!>

# Shift + Option + E  # to run code chunks
# Code completion: (Basic) Ctrl + Space  ;  (SmartType) Ctrl + Shift + Space

import os
import sys
import time
#import re
#import glob
#import pandas as pd

#print(sys.path)

# directory for raw fastq files
READDIR = '/scratch/pawsey1216/cliddicoat/ft2d_chn/1_meta_raw'

# set working directory
workDir = READDIR
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_index" % os.getcwd()
mkdir_p(job_directory)


def import_strings_from_file(filepath):
    """
    Imports a list of strings from a text file, with one string per line.
    Args:
        filepath (str): The path to the text file.
    Returns:
        list: A list of strings, where each string corresponds to a line
              from the file, with leading/trailing whitespace (including newlines) removed.
    """
    strings_list = []
    with open(filepath, 'r') as file:
        for line in file:
            strings_list.append(line.strip())
    return strings_list


file_path = '/scratch/pawsey1216/cliddicoat/ft2d_chn/1_meta_raw/forslund-t2d-chn-run-list.txt'
sra_runs = import_strings_from_file(file_path)
#print(sra_runs)



# iterate through creating jobs
n = len(sra_runs)
for i in range(n):
    #i=0
    job_file = os.path.join(job_directory, "%s_get_sra.sh" % sra_runs[i])
        
    #with open(job_file) as fh:
    fh = open(job_file, "w")
    fh.writelines("#!/bin/bash --login\n")
    fh.writelines("#SBATCH --account=pawsey1216\n")
    fh.writelines("#SBATCH --partition=work\n")
    #fh.writelines("#SBATCH --job-name=%s.job\n" % samples[i])
    #fh.writelines("#SBATCH --output=.out/%s.out\n" % samples[i])
    #fh.writelines("#SBATCH --error=.out/%s.err\n" % samples[i])
    fh.writelines("#SBATCH --ntasks=1\n")
    fh.writelines("#SBATCH --ntasks-per-node=1\n")
    fh.writelines("#SBATCH --cpus-per-task=24\n")
    fh.writelines("#SBATCH --time=24:00:00\n")
    #fh.writelines("#SBATCH --qos=normal\n")
    #fh.writelines("#SBATCH --mail-type=ALL\n")
    #fh.writelines("#SBATCH --mail-user=email_address\n")
    fh.writelines("export PATH=$PATH:/software/projects/pawsey1216/cliddicoat/sratoolkit.3.2.1-alma_linux64/bin\n")
    fh.writelines("fastq-dump --gzip --outdir %s --skip-technical --readids --read-filter pass --dumpbase --clip --split-3 %s\n" % (READDIR,sra_runs[i]))
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

## END
