# # # # # # # # # # # # # # # #
# # FastQC (without MultiQC) - for AMI 2025 metagenomics
# # # # # # # # # # # # # # # #
#$ sinfo -N -o "%N %C %e %m"   # get cluster usage info
#$ cd <check run folder!!>

# CHECK and ADJUST RANGE! 

import os
import sys
import time
import pandas as pd


print(sys.path)

# directory for raw fastq files
READDIR = '/scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw'

# directory for QC outputs
OUTDIR = '/scratch/pawsey1216/cliddicoat/ami_2025/1_qc'

# set working directory
workDir = OUTDIR
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_files" % os.getcwd()
mkdir_p(job_directory)

## read in manifest of sampleID and corresponding raw fastq R1/R2 files
# iterate through these with alternative filename patterns *_R1.fastq.gz, then *_R1_001.fastq.gz
#table = pd.read_csv('/scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw/combined_reads_dist_vs_nat.tsv.txt', sep='\t')
table = pd.read_csv('/scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw/combined_reads_001_dist_vs_nat.tsv.txt', sep='\t')

# extract info from pandas data table
r1_files = table["R1_filenames"]
r2_files = table["R2_filenames"]



# iterate through creating jobs for fastp QC of each sample # len(r1_files) = 285
n1 = len(r1_files)
n2 = len(r2_files)

#for i in range(n1):
#range(start,stop,step)
for i in range(1,950,50):
    #i=0
    job_file = os.path.join(job_directory, "slurm_fastqc_r1_%s.sh" % i)
    
    infile1 = os.path.join(READDIR,"%s" % r1_files[i] )
        
    #with open(job_file) as fh:
    fh = open(job_file, "w")
    fh.writelines("#!/bin/bash --login\n")
    fh.writelines("#SBATCH --account=pawsey1216\n")
    fh.writelines("#SBATCH --partition=work\n")
    fh.writelines("#SBATCH --ntasks=1\n")
    fh.writelines("#SBATCH --ntasks-per-node=1\n")
    fh.writelines("#SBATCH --cpus-per-task=16\n")
    fh.writelines("#SBATCH --time=2:00:00\n")
    
    fh.writelines("fastqc -o %s -t 16 %s\n" % (OUTDIR,infile1))
    
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)


#for i in range(n2):
##range(start,stop,step)
for i in range(1,950,50):
    #i=0
    job_file = os.path.join(job_directory, "slurm_fastqc_r2_%s.sh" % i)
    
    infile2 = os.path.join(READDIR,"%s" % r2_files[i] )
        
    #with open(job_file) as fh:
    fh = open(job_file, "w")
    fh.writelines("#!/bin/bash --login\n")
    fh.writelines("#SBATCH --account=pawsey1216\n")
    fh.writelines("#SBATCH --partition=work\n")
    fh.writelines("#SBATCH --ntasks=1\n")
    fh.writelines("#SBATCH --ntasks-per-node=1\n")
    fh.writelines("#SBATCH --cpus-per-task=16\n")
    fh.writelines("#SBATCH --time=2:00:00\n")
    
    fh.writelines("fastqc -o %s -t 16 %s\n" % (OUTDIR,infile2))
    
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

## END
