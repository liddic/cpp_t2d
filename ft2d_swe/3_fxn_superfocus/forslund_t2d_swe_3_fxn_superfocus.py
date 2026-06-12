# # # # # # # # # # # # # # # #
# # Run Superfocus fxn annotations from reads - for Forslund T2D Swedish metagenomics
# # # # # # # # # # # # # # # #
# Shift-Option-E  to run line or chunk
#$ sinfo -N -o "%N %C %e %m"   # get cluster usage info
#$ cd <check run folder!!>

import os
import sys
import time
import pandas as pd
import re


print(sys.path)


# qc reads and sample manifest table are here
READDIR = '/scratch/pawsey1216/cliddicoat/ft2d_swe/2b_clean_hostremoval'

# set working directory
workDir = '/scratch/pawsey1216/cliddicoat/ft2d_swe/3_fxn_superfocus'


# change to workDir
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


# set job directory for sbatch submissions
job_directory = "%s/job_files" % os.getcwd()
mkdir_p(job_directory)

# set temp directory
TEMP_DIR_START = '/scratch/pawsey1216/cliddicoat/temp'

## read in manifest of sampleID and corresponding raw fastq R1/R2 files
table = pd.read_csv('/scratch/pawsey1216/cliddicoat/ft2d_swe/1_meta_raw/combined_reads_ft2d_swe.tsv.txt', sep='\t')

# extract info from pandas data table
samples = table["sample"]
#r1_files = table["R1_filenames"]
#r2_files = table["R2_filenames"]

# iterate through creating jobs
n = len(samples)

for i in range(n):
    #i=0
    job_file = os.path.join(job_directory, "submission_superfocus_%s.sh" % list(samples)[i])

    #qc_read_r1 = os.path.join(READDIR, "%s_R1.good.fastq" % list(samples)[i])
    qc_read_r1 = os.path.join(READDIR,"%s_non_host.1.fastq" % list(samples)[i])

    OUTDIR = "%s/superfocus_out_%s" % (os.getcwd(),list(samples)[i])
    #mkdir_p(OUTDIR)  # this output directory will be made by superfocus

    path_to_db_dir = "/software/projects/pawsey1216/cliddicoat/conda_envs/_metg_fxn_superfocus/lib/python3.8/site-packages/superfocus-0.0.0-py3.8.egg/superfocus_app"

    temp_directory = "%s/%s" % (TEMP_DIR_START, list(samples)[i])
    mkdir_p(temp_directory)

    # write and launch submission scripts
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
    #fh.writelines("#SBATCH --nodelist=hpc-node0XX,hpc-node0XX,hpc-node0XX,hpc-node0XX,hpc-node0XX,hpc-node0XX\n")
    #fh.writelines("#SBATCH --nodelist=hpc-node019,hpc-node020,hpc-node021\n")
    #fh.writelines("#SBATCH --qos=normal\n")
    #fh.writelines("#SBATCH --mail-type=ALL\n")
    #fh.writelines("#SBATCH --mail-user=email_address\n")

    # set the TMPDIR environment variable
    fh.writelines("export TMPDIR=%s\n" % (temp_directory))
    fh.writelines("export PATH=$PATH:/software/projects/pawsey1216/cliddicoat/conda_envs/_metg_fxn_superfocus/bin\n")

    ### Run Superfocus - just on R1 reads

    # superfocus -q <path-to-fasta/q-files-or-directory> -dir <path-to-output-directory> -a <aligner> -db DATABASE
    fh.writelines("superfocus -t 32 -q %s -dir %s -a diamond -db DB_100 --alternate_directory %s\n" % (qc_read_r1,OUTDIR,path_to_db_dir))
    # -b --alternate_directory # Alternate directory for your databases

    fh.close()

    time.sleep(2)  # Wait between submitting jobs

    os.system("sbatch %s" % job_file) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

## END