#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw
#find ./ -type f -name "*.fastq.gz"
find ./ -type f -name "*.fastq.gz" -exec md5sum "{}" + > ami_dist_vs_nat_md5sum_checklist_gz.txt
