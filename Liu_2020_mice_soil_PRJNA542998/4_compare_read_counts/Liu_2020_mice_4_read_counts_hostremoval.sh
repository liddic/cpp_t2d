#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/4_compare_read_counts

awk '{print $1}' /scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/1_meta_raw/combined_reads_Liu_2020_mice_soil_PRJNA542998.tsv.txt > samples3.txt

## remove header
awk 'NR>1' samples3.txt > samples_noheader3.txt

for f in `cat samples_noheader3.txt`; do line_count=$(wc -l < /scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/2b_clean_hostremoval/"$f"_non_host.1.fastq) && result=$(($line_count / 4)) && echo "$f $result"; done > nonhost_r1_sequence_counts.txt
