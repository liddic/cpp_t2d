#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/ft2d_chn/2c_rarefy_even_seq_min
export PATH=$PATH:/software/projects/pawsey1216/cliddicoat/conda_envs/_seqtk/bin

awk '{print $1}' /scratch/pawsey1216/cliddicoat/ft2d_chn/2c_rarefy_even_seq_min/keep_t2d_chn_list_min.txt > samples_noheader.txt

for f in `cat samples_noheader.txt`; do

# seqtk sample -s100 input.fastq 5000 > output.fastq
# CHN Min
# 758416
# total n = 82

seqtk sample -s100 /scratch/pawsey1216/cliddicoat/ft2d_chn/2b_clean_hostremoval/"$f"_non_host.1.fastq 758416 > /scratch/pawsey1216/cliddicoat/ft2d_chn/2c_rarefy_even_seq_min/"$f"_non_host_rarefy_even.1.fastq

done
