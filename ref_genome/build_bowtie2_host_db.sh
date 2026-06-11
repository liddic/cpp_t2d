#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00

export PATH=$PATH:/software/projects/pawsey1216/cliddicoat/conda_envs/_metg_bowtie2/bin

cd /scratch/pawsey1216/cliddicoat/ref_genome

# bowtie2-build host_genome.fna host_DB

bowtie2-build GCF_000001405.40_GRCh38.p14_genomic.fna.gz host_db
