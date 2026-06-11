#!/bin/bash --login
#SBATCH --account=pawsey1216
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

cd /scratch/pawsey1216/cliddicoat/sunbad_resto/cleaned

/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679658.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679659.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679660.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679661.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679662.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679663.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679664.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679665.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679666.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679667.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679668.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679669.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679670.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679671.3 --file 299.1
/software/projects/pawsey1216/cliddicoat/miniconda3/bin/mg-download.py --metagenome mgm4679672.3 --file 299.1
