#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/ami_2025/4_read_counts

find /scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw -type f -name "*_R1.fastq.gz" -print | sort > read1_files.txt
for f in `cat read1_files.txt`; do line_count=$(zcat "$f" | wc -l) && result=$(($line_count / 4)) && echo "$f $result"; done > r1_meta_raw_seq_counts_zcat.txt

find /scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw -type f -name "*_R2.fastq.gz" -print | sort > read2_files.txt
for f in `cat read2_files.txt`; do line_count=$(zcat "$f" | wc -l) && result=$(($line_count / 4)) && echo "$f $result"; done > r2_meta_raw_seq_counts_zcat.txt

# 001 files

find /scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw -type f -name "*_R1_001.fastq.gz" -print | sort > read1_001_files.txt
for f in `cat read1_001_files.txt`; do line_count=$(zcat "$f" | wc -l) && result=$(($line_count / 4)) && echo "$f $result"; done > r1_meta_raw_001_seq_counts_zcat.txt

find /scratch/pawsey1216/cliddicoat/ami_2025/1_meta_raw -type f -name "*_R2_001.fastq.gz" -print | sort > read2_001_files.txt
for f in `cat read2_001_files.txt`; do line_count=$(zcat "$f" | wc -l) && result=$(($line_count / 4)) && echo "$f $result"; done > r2_meta_raw_001_seq_counts_zcat.txt


# END
