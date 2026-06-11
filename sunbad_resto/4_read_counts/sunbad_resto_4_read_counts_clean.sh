#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/sunbad_resto/4_read_counts

awk '{print $1}' /scratch/pawsey1216/cliddicoat/sunbad_resto/2_clean_mgrast/sunbad-resto-mgp16379_metadata.tsv.txt > samples.txt

## remove header
awk 'NR>1' samples.txt > samples_noheader.txt

# fasta fna files so divide line count by 2

for f in `cat samples_noheader.txt`; do line_count=$(wc -l < /scratch/pawsey1216/cliddicoat/sunbad_resto/2_clean_mgrast/"$f".299.screen.passed.fna) && result=$(($line_count / 2)) && echo "$f $result"; done > clean_r1_sequence_counts.txt
