# # # # # # # # # # # # # # # #
# # Download AMI fastq files - Soils x Metagenomics Disturbed vs Natural comparison dataset
# # # # # # # # # # # # # # # #
#$ sinfo -N -o "%N %C %e %m"   # get cluster usage info
#$ cd <check run folder!!>

# Shift + Option + E  # to run code chunks
# Code completion: (Basic) Ctrl + Space  ;  (SmartType) Ctrl + Shift + Space

import os
import sys
import time
import requests
#import hashlib   # for md5sum check downloads

# set working directory
READDIR = '/scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw'
os.chdir(READDIR)

## read in list of urls
with open( "ami_soil_metagenomics_url_list_disturbed_vs_natural_dataset.txt" , "r" ) as f:
    urls = f.readlines()


with open( "/home/cliddicoat/CKAN_API_KEY__CL.txt" , "r" ) as f:
    CKAN_API_KEY = f.read()


def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)


job_directory = "%s/job_index" % os.getcwd()
mkdir_p(job_directory)


# iterate through creating jobs
n = len(urls)
#for i in range(2):
for i in range(n):
#for i in range(1356,1929):

    #i=0
    
    thisDownload = urls[i]
    thisDownload = thisDownload.replace("\n", "")
    thisFileName = os.path.basename(thisDownload)
    thisNo = i
    
    python_file = os.path.join(job_directory, "get_ami_%s.py" % thisNo)
    fh = open(python_file, "w")
    fh.writelines("import os\n")
    fh.writelines("import sys\n")
    fh.writelines("import time\n")
    fh.writelines("import requests\n")
    fh.writelines("READDIR = '/scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw'\n")
    fh.writelines("os.chdir(READDIR)\n")
    fh.writelines("with open( '/home/cliddicoat/CKAN_API_KEY__CL.txt' , 'r' ) as f:\n")
    fh.writelines("    CKAN_API_KEY = f.read()\n")
    fh.writelines("\n")
    fh.writelines("\n")
    fh.writelines("HEADERS = {'Authorization':CKAN_API_KEY ,'Destination':READDIR }\n")
    fh.writelines("\n")
    fh.writelines("res = requests.request('GET', '%s' , headers = HEADERS , allow_redirects=True)\n" % thisDownload)
    fh.writelines("\n")
    fh.writelines("open( os.path.join(READDIR, 'x%s') ,'wb').write(res.content)\n" % thisFileName)
    fh.writelines("\n")
    fh.close()
    time.sleep(2)
    
    job_file = os.path.join(job_directory, "ami_%s.sh" % thisNo)
        
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
    fh.writelines("#SBATCH --cpus-per-task=16\n")
    fh.writelines("#SBATCH --time=24:00:00\n")
    #fh.writelines("#SBATCH --qos=normal\n")
    #fh.writelines("#SBATCH --mail-type=ALL\n")
    #fh.writelines("#SBATCH --mail-user=email_address\n")
    fh.writelines("python %s\n" % python_file)
    fh.close()
    time.sleep(2)
    os.system("sbatch %s" % job_file)

## END