#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00

module load r/4.4.1

#export R_LIBS_USER=/software/projects/<project-id>/<user-name>/setonix/<DATE-TAG>/r/%v
#export R_LIBS_USER=/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1

cd /scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th

# Execute the R script

Rscript Rcode_cpp3d_onPawsey_t2d_swe_rarefied_15th.R
