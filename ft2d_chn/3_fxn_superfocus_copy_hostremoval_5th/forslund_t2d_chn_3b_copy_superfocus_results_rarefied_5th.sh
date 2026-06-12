#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/ft2d_chn/3_fxn_superfocus_copy_hostremoval_5th

awk '{print $1}' /scratch/pawsey1216/cliddicoat/ft2d_chn/2c_rarefy_even_seq_5th/keep_t2d_chn_list_5th.txt > samples_noheader.txt

## remove header
#awk 'NR>1' samples.txt > samples_noheader.txt

for f in `cat samples_noheader.txt`; do

mkdir -p superfocus_out_"$f"
cp /scratch/pawsey1216/cliddicoat/ft2d_chn/3_fxn_superfocus_5th/superfocus_out_"$f"/output_all_levels_and_function.xls /scratch/pawsey1216/cliddicoat/ft2d_chn/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_"$f"/output_all_levels_and_function.xls

done
