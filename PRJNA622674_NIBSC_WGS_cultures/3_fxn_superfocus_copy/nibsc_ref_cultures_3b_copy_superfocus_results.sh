#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/3_fxn_superfocus_copy

awk '{print $1}' /scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/1_meta_raw/combined_reads_PRJNA622674_NIBSC_WGS_cultures.tsv.txt > samples.txt

# remove header
awk 'NR>1' samples.txt > samples_noheader.txt

for f in `cat samples_noheader.txt`; do

mkdir -p superfocus_out_"$f"
cp /scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/3_fxn_superfocus/superfocus_out_"$f"/output_all_levels_and_function.xls /scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/3_fxn_superfocus_copy/superfocus_out_"$f"/output_all_levels_and_function.xls

done
