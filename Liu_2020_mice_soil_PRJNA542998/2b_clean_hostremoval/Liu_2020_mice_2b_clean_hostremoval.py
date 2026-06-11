# # # # # # # # # # # # # # # #
# # Bowtie2 - human host genome sequence removal - for Liu_2020_mice_soil_PRJNA542998  metagenomics
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
import pandas as pd


print(sys.path)

# directory for input fastq files
READDIR = '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/2_clean_fastp'

# directory for non-host outputs
OUTDIR = '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/2b_clean_hostremoval'

# set working directory
workDir = OUTDIR
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_index" % os.getcwd()
mkdir_p(job_directory)

## read in manifest of sampleID and corresponding raw fastq R1/R2 files
table = pd.read_csv('/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/1_meta_raw/combined_reads_Liu_2020_mice_soil_PRJNA542998.tsv.txt', sep='\t')

# extract info from pandas data table
samples = table["sample"]
#r1_files = table["R1_filenames"]
#r2_files = table["R2_filenames"]

# iterate through creating jobs for fastp QC of each sample
n = len(samples)
for i in range(n):
    #i=0
    job_file = os.path.join(job_directory, "%s.sh" % samples[i])
    
    trimmed_r1 = os.path.join(READDIR,"%s_R1.good.fastq" % samples[i] )
    trimmed_r2 = os.path.join(READDIR,"%s_R2.good.fastq" % samples[i] )
    
    output_file = os.path.join(OUTDIR,"%s_non_host.fastq" % samples[i] )
    
    threads = 32
    #db_path = "/scratch/pawsey1216/cliddicoat/ref_genome/host_db"
    db_path = "/scratch/pawsey1216/cliddicoat/ref_genome_MOUSE/host_mouse_db"
        
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
    fh.writelines("#SBATCH --cpus-per-task=32\n")
    fh.writelines("#SBATCH --time=24:00:00\n")
    #fh.writelines("#SBATCH --qos=normal\n")
    #fh.writelines("#SBATCH --mail-type=ALL\n")
    #fh.writelines("#SBATCH --mail-user=email_address\n")
    fh.writelines("export PATH=$PATH:/software/projects/pawsey1216/cliddicoat/conda_envs/_metg_bowtie2/bin\n")
    fh.writelines("bowtie2 -x %s -1 %s -2 %s -p %s --un-conc %s \n" % (db_path,trimmed_r1,trimmed_r2,threads,output_file))
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

## END
