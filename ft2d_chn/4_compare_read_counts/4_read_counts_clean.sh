#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/ft2d_chn/4_compare_read_counts

awk '{print $1}' /scratch/pawsey1216/cliddicoat/ft2d_chn/1_meta_raw/forslund-t2d-chn-run-list.txt > samples_noheader.txt

for f in `cat samples_noheader.txt`; do line_count=$(wc -l < /scratch/pawsey1216/cliddicoat/ft2d_chn/2_clean_fastp/"$f"_R1.good.fastq) && result=$(($line_count / 4)) && echo "$f $result"; done > clean_r1_sequence_counts.txt
