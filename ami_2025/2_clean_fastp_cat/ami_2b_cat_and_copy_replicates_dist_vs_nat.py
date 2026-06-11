# # # # # # # # # # # # # # # #
# # concatenate and copy replicate cleaned fastq files
# # # # # # # # # # # # # # # #

import os
import sys
import time
import re

print(sys.path)

# set working directory
workDir = '/scratch/pawsey1216/cliddicoat/ami_2025/2_clean_fastp_cat'
DIR_READS= '/scratch/pawsey1216/cliddicoat/ami_2025/2_clean_fastp'
OUTDIR = '/scratch/pawsey1216/cliddicoat/ami_2025/2_clean_fastp_cat'


# change to workDir
os.chdir(workDir)
print("Current Working Directory: ",os.getcwd())

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_cat_copy" % os.getcwd()
mkdir_p(job_directory)

## read in list of samples
with open( "ami_dist_vs_nat_unique_samps.txt" , "r" ) as f:
    samples = [line.strip() for line in f.readlines()]


## identify the read files
qcreadsFileNames = [ name for name in os.listdir(DIR_READS) if os.path.isfile(os.path.join(DIR_READS, name)) ]

# only keep .fastq files
reads = [i for i in qcreadsFileNames if i.endswith('.fastq')]

# convert to full paths
reads_fullpath = []
for element in reads:
    reads_fullpath.append(os.path.join(DIR_READS, element))


n = len(samples)

for i in range(n):
    job_file = os.path.join(job_directory, "%s.sh" % samples[i])
    
    # need to define: input_r1, input_r2, input_s
    reg = re.compile(r'%s' % samples[i])  # Compile the regex
    samp_reads = list(filter(reg.search, reads_fullpath))  # Create iterator using filter, cast to list

    reg1 = re.compile(r'_R1.good')  # Compile the regex
    input_r1 = list(filter(reg1.search, samp_reads))
    input_r1 = " ".join(input_r1)

    reg2 = re.compile(r'_R2.good')  # Compile the regex
    input_r2 = list(filter(reg2.search, samp_reads))
    input_r2 = " ".join(input_r2)

    regS1 = re.compile(r'_R1.single')  # Compile the regex
    input_s1 = list(filter(regS1.search, samp_reads))
    input_s1 = " ".join(input_s1)
    
    regS2 = re.compile(r'_R2.single')  # Compile the regex
    input_s2 = list(filter(regS2.search, samp_reads))
    input_s2 = " ".join(input_s2)
        
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
    fh.writelines("#SBATCH --cpus-per-task=16\n")
    fh.writelines("#SBATCH --time=24:00:00\n")
    #fh.writelines("#SBATCH --qos=normal\n")
    #fh.writelines("#SBATCH --mail-type=ALL\n")
    #fh.writelines("#SBATCH --mail-user=email_address\n")
    fh.writelines("cat %s > %s/%s_R1.good.fastq\n" % (input_r1, OUTDIR, samples[i] ))
    fh.writelines("cat %s > %s/%s_R2.good.fastq\n" % (input_r2, OUTDIR, samples[i] ))
    fh.writelines("cat %s > %s/%s_R1.single.fastq\n" % (input_s1, OUTDIR, samples[i] ))
    fh.writelines("cat %s > %s/%s_R2.single.fastq\n" % (input_s2, OUTDIR, samples[i] ))
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

# END