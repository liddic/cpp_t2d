#########################
#
# cpp_t2d: Examining shared trends in compound-associated functional capacities 
# of degraded ecosystem soil microbiomes and type 2 diabetes gut microbiomes
#
# Using compound processing potential (CPP): https://github.com/liddic/cpp_t2d
# Craig Liddicoat | Flinders University, South Australia 
#
# PART 3 - R code for supplementary data analyses
#
#########################

# record library and version info
.libPaths() # "/Library/Frameworks/R.framework/Versions/4.2/Resources/library"

R.Version()
# "R version 4.2.2 (2022-10-31)"
citation()
# R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
# https://www.R-project.org/.

library(readxl); packageVersion("readxl") # '1.4.2'
library(plyr); packageVersion("plyr") # '1.8.8'
library(dplyr); packageVersion("dplyr") # ‘1.1.2’
library(vegan);packageVersion("vegan") # '2.6.4'
library(phyloseq); packageVersion("phyloseq") # ‘1.44.0’
library(ggplot2); packageVersion("ggplot2") # ‘3.5.1’
library(grid); packageVersion("grid") #  '4.2.2'
library(reshape2); packageVersion("reshape2") # '1.4.4'
library(tidyr); packageVersion("tidyr") # ‘1.3.0’
library(corrr); packageVersion("corrr") # '0.4.4'
library(ggforce); packageVersion("ggforce") # '0.4.1'
library(ggrepel); packageVersion("ggrepel") # '0.9.2'
library(stringdist); packageVersion("stringdist") # ‘0.9.10’
library(stringr); packageVersion("stringr") # ‘1.5.0’
library(doParallel); packageVersion("doParallel") # '1.0.17'
library(RColorBrewer); packageVersion("RColorBrewer") # '1.1.3'
library(ggpp); packageVersion("ggpp") # ‘0.5.0’
library(MASS); packageVersion("MASS") # ‘7.3.58.1’
library(ggsignif); packageVersion("ggsignif") # '0.6.4'
library(moments); packageVersion("moments") # ‘0.14.1’
library(grDevices); packageVersion("grDevices") #  '4.2.2'
library(ggbiplot); packageVersion("ggbiplot") #  ‘0.55’
library(viridis); packageVersion("viridis") #  ‘0.6.3’
library(FSA); packageVersion("FSA") # '0.9.3'
library(rcompanion); packageVersion("rcompanion") # '2.4.18'
library(fields); packageVersion("fields") # ‘14.1’
library(car); packageVersion("car") # ‘3.1.1’
library(multcompView); packageVersion("multcompView") # ‘0.1.8’
library(gtools); packageVersion("gtools") # ‘3.9.4’
library(igraph); packageVersion("igraph") #  '1.4.2'
library(pheatmap); packageVersion("pheatmap") # '1.0.12'
library(RColorBrewer); packageVersion("RColorBrewer") # ‘1.1.3’
library(VennDiagram); packageVersion("VennDiagram") # ‘1.7.3’


#########################
## CAUTION!! save.image("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v8h.RData")
##      load("/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R/cpp3d-indiv-resto-vs-t2d-WORKSPACE-v8h.RData")
#########################

workdir <- "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"
setwd(workdir)
getwd()

par.default <- par()

##########################
##########################
##########################
##########################

## Supplementary validation datasets

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - Liu_2020_mice
#    prepare for data download from SRA
#-------------------------

# Liu_2020_mice_soil_PRJNA542998
# 
# Exposure to soil environments during earlier life stages is distinguishable in the gut microbiome of adult mice
#
# Wenjun Liu , Zheng Sun , Chen Ma , Jiachao Zhang , ChenChen Ma , Yinqi
# Zhao , Hong Wei , Shi Huang & Heping Zhang (2021) Exposure to soil environments during earlier
# life stages is distinguishable in the gut microbiome of adult mice, Gut Microbes, 13:1, 1-13, DOI:
#   10.1080/19490976.2020.1830699
# 
# Totally, 90 germ-free mice were randomized
# into three groups (30 per group) and each
# group was raised under an assigned-simulated
# environment for 60 days, using different natural
# soil samples collected from steppe, forest or desert
# habitats.
# 
# Mice_soil_exposure
# desert / steppe / forest
# Compare steppe versus forest 
# 
# From Biosample - extract 
# 
# F_first_sampling (forest)
# G_first_sampling (grassland)
# D_first_sampling (desert)


sradat <- read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/SraRunTable_Liu_2020_mice_soil_PRJNA542998_LinkedData.xlsx", sheet = 1, range = "A1:AH202")
sradat <- as.data.frame(sradat)
str(sradat)
# 'data.frame':	201 obs. of  34 variables:
# $ Run                           : chr  "SRR9276250" "SRR9276308" "SRR9276307" "SRR9276316" ...
# $ Assay Type                    : chr  "WGS" "WGS" "WGS" "WGS" ...
# $ AvgSpotLen                    : num  252 252 252 252 568 252 252 252 252 252 ...
# $ Bases                         : num  8.78e+09 8.59e+09 7.81e+09 8.76e+09 8.53e+09 ...
# $ BioProject                    : chr  "PRJNA542998" "PRJNA542998" "PRJNA542998" "PRJNA542998" ...
# $ BioSample                     : chr  "SAMN11664431" "SAMN11664430" "SAMN11664429" "SAMN11664428" ...
# $ BioSampleModel                : chr  "Metagenome or environmental" "Metagenome or environmental" "Metagenome or environmental" "Metagenome or environmental" ...
# $ Bytes                         : num  3.35e+09 3.28e+09 2.97e+09 3.33e+09 5.36e+09 ...
# $ Center Name                   : chr  "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" ...
# $ Collection_Date               : POSIXct, format: "2014-12-01" "2014-12-01" "2014-12-01" "2014-12-01" ...
# $ Consent                       : chr  "public" "public" "public" "public" ...
# $ DATASTORE filetype            : chr  "fastq,run.zq,sra" "run.zq,fastq,sra" "fastq,sra,run.zq" "sra,run.zq,fastq" ...
# $ DATASTORE provider            : chr  "ncbi,s3,gs" "s3,gs,ncbi" "ncbi,s3,gs" "gs,s3,ncbi" ...
# $ DATASTORE region              : chr  "gs.us-east1,ncbi.public,s3.us-east-1" "ncbi.public,gs.us-east1,s3.us-east-1" "gs.us-east1,s3.us-east-1,ncbi.public" "gs.us-east1,ncbi.public,s3.us-east-1" ...
# $ Experiment                    : chr  "SRX6046069" "SRX6046011" "SRX6046012" "SRX6046003" ...
# $ geo_loc_name_country          : chr  "China" "China" "China" "China" ...
# $ geo_loc_name_country_continent: chr  "Asia" "Asia" "Asia" "Asia" ...
# $ geo_loc_name                  : chr  "China" "China" "China" "China" ...
# $ HOST                          : chr  "mice" "mice" "mice" "mice" ...
# $ Instrument                    : chr  "Illumina HiSeq 2500" "Illumina HiSeq 2500" "Illumina HiSeq 2500" "Illumina HiSeq 2500" ...
# $ isolation_source              : chr  "mice fecals" "mice fecals" "mice fecals" "mice fecals" ...
# $ lat_lon                       : chr  "29.35 N 106.33 E" "29.35 N 106.33 E" "29.35 N 106.33 E" "29.35 N 106.33 E" ...
# $ Library Name                  : chr  "D13440" "D13439" "D13438" "D13437" ...
# $ LibraryLayout                 : chr  "PAIRED" "PAIRED" "PAIRED" "PAIRED" ...
# $ LibrarySelection              : chr  "DNase" "DNase" "DNase" "DNase" ...
# $ LibrarySource                 : chr  "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" ...
# $ Organism                      : chr  "feces metagenome" "feces metagenome" "feces metagenome" "feces metagenome" ...
# $ Platform                      : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
# $ ReleaseDate                   : chr  "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" ...
# $ create_date                   : chr  "2019-06-11T22:12:00Z" "2019-06-11T22:31:00Z" "2019-06-11T22:19:00Z" "2019-06-11T22:26:00Z" ...
# $ version                       : num  1 1 1 1 1 1 1 1 1 1 ...
# $ Sample Name                   : chr  "D13440" "D13439" "D13438" "D13437" ...
# $ SRA Study                     : chr  "SRP201145" "SRP201145" "SRP201145" "SRP201145" ...
# $ Group                         : chr  "F2F_second_sampling" "F2F_second_sampling" "F2F_second_sampling" "F2F_second_sampling" ...

length(unique(sradat$`Sample Name`)) # 201

unique(sradat$BioProject) # "PRJNA542998"

table(sradat$Group)
# D_first_sampling   D_taking_environment_samples            D2D_second_sampling D2D_taking_environment_samples 
# 30                              4                             10                              1 
# D2F_second_sampling D2F_taking_environment_samples            D2G_second_sampling D2G_taking_environment_samples 
# 10                              1                             10                              1 
# F_first_sampling   F_taking_environment_samples            F2D_second_sampling F2D_taking_environment_samples 
# 30                              4                             10                              1 
# F2F_second_sampling F2F_taking_environment_samples            F2G_second_sampling F2G_taking_environment_samples 
# 10                              1                             10                              1 
# G_first_sampling   G_taking_environment_samples            G2D_second_sampling G2D_taking_environment_samples 
# 30                              4                             10                              1 
# G2F_second_sampling G2F_taking_environment_samples            G2G_second_sampling G2G_taking_environment_samples 
# 10                              1                             10                              1 

# only consider 1W (=T2D), versus 2W (T2D with 3 months high fibre/ WTP diet intervention)

sel <- which(sradat$Group %in% c("D_first_sampling", "G_first_sampling" , "F_first_sampling",
                                 "D_taking_environment_samples", "G_taking_environment_samples" , "F_taking_environment_samples" )) # 102
30*3 + 4*3 # 102

sradat.select <- sradat[sel, ]

dim(sradat.select) # 102  34

table(sradat.select$Group)
# D_first_sampling D_taking_environment_samples             F_first_sampling F_taking_environment_samples             G_first_sampling 
# 30                            4                           30                            4                           30 
# G_taking_environment_samples 
# 4 


temp <- sradat.select

saveRDS(object = sradat.select, file = "sradat.select.Liu_2020_mice_soil_PRJNA542998.rds")


liu_runs <- sradat.select$Run
length(liu_runs) # 102

# file for SRA download

writeLines(liu_runs, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/Liu_2020_mice-run-list.txt")

#write.csv(x = paste0("fastq-dump --outdir /scratch/user/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw --skip-technical --readids --read-filter pass --dumpbase --clip --split-3 ",sort(chn_runs)), file = "forslund-t2d-CHN-testset-sra-runs-download.txt", quote = FALSE, row.names = FALSE)
#write.csv(x = paste0("fastq-dump --outdir /cluster/jobs/lidd0026/forslund-t2d-chn/ft2d_1_meta_raw --skip-technical --readids --read-filter pass --dumpbase --clip --split-3 ",sort(chn_runs)), file = "forslund-t2d-CHN-testset-sra-runs-download-cluster.txt", quote = FALSE, row.names = FALSE)


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - Liu_2020_mice - w/ Host-removal - compare no reads w/ vs w/out host-removal
#-------------------------

r1.clean <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/4_compare_read_counts/clean_r1_sequence_counts.txt", header = FALSE, sep = " ")
head(r1.clean)
names(r1.clean) <- c("sample","n_reads_cleaned")

r1.nonhost <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/4_compare_read_counts/nonhost_r1_sequence_counts.txt", header = FALSE, sep = " ")
head(r1.nonhost)
names(r1.nonhost) <- c("sample","n_reads_nonhost")

identical(r1.clean$sample, r1.nonhost$sample) # TRUE

r1.mice <- r1.clean
r1.mice$n_reads_nonhost <- r1.nonhost$n_reads_nonhost
r1.mice$percent_host_removed <- 100*(r1.mice$n_reads_cleaned - r1.mice$n_reads_nonhost)/r1.mice$n_reads_cleaned
hist(r1.mice$percent_host_removed);summary(r1.mice$percent_host_removed)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00091  0.20162  0.47056  1.82501  1.27338 43.64327 

r1.mice$n_reads_host_removed <- r1.mice$n_reads_cleaned - r1.mice$n_reads_nonhost

median(r1.mice$percent_host_removed) # 0.4705625

# For the Liu_2020_mice_soil_PRJNA542998 samples contained a median (and interquartile range, IQR): 0.47056 % (0.20162 - 1.27338 %) of reads that were classified as mouse sequences.

hist(r1.mice$n_reads_cleaned);summary(r1.mice$n_reads_cleaned)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 7976674 21167172 24302852 26480025 31807341 55420686 

hist(r1.mice$n_reads_nonhost);summary(r1.mice$n_reads_nonhost)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 7784072 20467986 24130459 26051635 31266072 55411845 

# For the Liu_2020_mice_soil_PRJNA542998 data we used median 24130459 (IQR: 204679868 - 31266072) cleaned non-host sequences

hist(r1.mice$n_reads_host_removed);summary(r1.mice$n_reads_host_removed)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  382   38632  131801  428390  339397 9724941 

# median 131801 (IQR 38632 - 339397)

#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host-removal - read in superfocus - fxn potential outputs
#-------------------------

#saveRDS(object = sradat.select, file = "sradat.select.Liu_2020_mice_soil_PRJNA542998.rds")

sradat.select <- readRDS("sradat.select.Liu_2020_mice_soil_PRJNA542998.rds")

sampid <- sradat.select$Run # 102

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276169"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276170"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276171"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276174"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276175"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276176"


names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR9276169 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276169" 
# SRR9276170 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276170" 
# SRR9276171 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276171" 
# SRR9276174 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276174" 
# SRR9276175 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276175" 
# SRR9276176 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/3_fxn_superfocus_copy_hostremoval/superfocus_out_SRR9276176" 


length(results_dirs) # 101

sel <- which(names(results_dirs) %in% sampid) # qty 101
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 101

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # FALSE
length(results_dirs) # 101
length(sampid) # 102 - one sample did not have matching R1/R2 read in SRA repository

# reset sampid to remove missing sample
sampid <- names(results_dirs)
identical(sampid, names(results_dirs)) # TRUE


# In this data one Run corresponds to a single Sample_ID !!!

# collate results into a long-format table

sfx.long <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,percent_abun=NA)

for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  # [1] "Subsystem.Level.1"                                                                         
  # [2] "Subsystem.Level.2"                                                                         
  # [3] "Subsystem.Level.3"                                                                         
  # [4] "Function"                                                                                  
  # [5] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq"  
  # [6] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq.."
  
  
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [7] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [8] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [9] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [10] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  # [11] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [12] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.fastq..$", x = names(tab)) # for single (unpaired) reads
  if (length(sel.col.percent)>1) {
    tab$percent_abun <- apply(X = tab[ ,sel.col.percent], MARGIN = 1, FUN = mean )
  } else {
    tab$percent_abun <- tab[ ,sel.col.percent]
  }
  
  # sum(tab$percent_abun) # 100
  # mean(tab$percent_abun) # 0.004338583
  
  names(sfx.long) # "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"    "percent_abun"
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # ...
  # [13] "sampid"
  # [14] "percent_abun"
  
  tab <- tab[ ,c("sampid","Subsystem.Level.1","Subsystem.Level.2","Subsystem.Level.3","Function","percent_abun")]
  names(tab) <- names(sfx.long)
  
  sfx.long <- rbind(sfx.long,tab)
  
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}


head(sfx.long)
# remove empty 1st row
sfx.long <- sfx.long[-1, ]
dim(sfx.long) # 1232593       6
head(sfx.long)
#     sampleID                   subsys_L1 subsys_L2           subsys_L3                                                                                  fxn percent_abun
# 2 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase                                                2-methylaconitate_cis-trans_isomerase 6.361540e-05
# 3 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase        2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117) 2.226539e-05
# 4 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase                                              4-hydroxyproline_epimerase_(EC_5.1.1.8) 3.602222e-03
# 5 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase                                                        Alanine_racemase_(EC_5.1.1.1) 6.997694e-04
# 6 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase                                        Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic 4.771155e-06
# 7 SRR9276169 Amino Acids and Derivatives         - Amino acid racemase Alanine_racemase_(EC_5.1.1.1)_#_present_in_exosporium,_involved_in_spore_germination 2.775222e-04

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 31719   102

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 31719
length(unique(full_fxn_names)) # 31719

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_isomerase" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-oxalomesaconate_tautomerase_(EC_5.3.2.8)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR9276169"   "SRR9276170"   "SRR9276171"   "SRR9276174"   "SRR9276175"   "SRR9276176"   "SRR9276177"   "SRR9276178"   "SRR9276179"   "SRR9276180"  
# [12] "SRR9276181"   "SRR9276182"   "SRR9276183"   "SRR9276184"   "SRR9276185"   "SRR9276186"   "SRR9276188"   "SRR9276189"   "SRR9276190"   "SRR9276191"   "SRR9276192"  
# [23] "SRR9276193"   "SRR9276196"   "SRR9276197"   "SRR9276224"   "SRR9276225"   "SRR9276226"   "SRR9276227"   "SRR9276228"   "SRR9276229"   "SRR9276230"   "SRR9276231"  
# [34] "SRR9276232"   "SRR9276233"   "SRR9276243"   "SRR9276244"   "SRR9276245"   "SRR9276246"   "SRR9276247"   "SRR9276248"   "SRR9276249"   "SRR9276251"   "SRR9276252"  
# [45] "SRR9276253"   "SRR9276254"   "SRR9276255"   "SRR9276256"   "SRR9276257"   "SRR9276258"   "SRR9276259"   "SRR9276260"   "SRR9276269"   "SRR9276270"   "SRR9276271"  
# [56] "SRR9276272"   "SRR9276275"   "SRR9276277"   "SRR9276278"   "SRR9276279"   "SRR9276280"   "SRR9276281"   "SRR9276294"   "SRR9276303"   "SRR9276304"   "SRR9276305"  
# [67] "SRR9276306"   "SRR9276317"   "SRR9276318"   "SRR9276319"   "SRR9276320"   "SRR9276321"   "SRR9276322"   "SRR9276323"   "SRR9276324"   "SRR9276325"   "SRR9276326"  
# [78] "SRR9276335"   "SRR9276336"   "SRR9276337"   "SRR9276338"   "SRR9276339"   "SRR9276340"   "SRR9276341"   "SRR9276342"   "SRR9276343"   "SRR9276344"   "SRR9276345"  
# [89] "SRR9276346"   "SRR9276347"   "SRR9276356"   "SRR9276359"   "SRR9276360"   "SRR9276361"   "SRR9276362"   "SRR9276363"   "SRR9276364"   "SRR9276365"   "SRR9276366"  
# [100] "SRR9276367"   "SRR9276368"   "SRR9276369"  

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR9276169" "SRR9276170" "SRR9276171" "SRR9276174" "SRR9276175" "SRR9276176"

length(sampid) # 101

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # TRUE
#identical(sort(names(sfx.wide)), sort(as.character(sampid))) #

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 31719     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 194
length(unique(tax.fxn$subsys_L3)) # 1265
length(unique(tax.fxn$fxn)) # 15824


# # # #

## gather Function count data??
sfx.long.count <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,count_abun=NA)
length(sampid) # 101
for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  tab$sampid <- this_samp
  names(tab)
  tab <- tab[,c(7,1,2,3,4,5)] # this time capture 'count' data
  names(tab) <- names(sfx.long.count)
  sfx.long.count <- rbind(sfx.long.count,tab)
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}
head(sfx.long.count)
# remove empty 1st row
sfx.long.count <- sfx.long.count[-1, ]
sum(sfx.long.count$count_abun) # 261705854 = 261,705,854
sfx.long.count$full_fxn_tax <- paste0(sfx.long.count$subsys_L1,"___", sfx.long.count$subsys_L2,"___", sfx.long.count$subsys_L3,"___", sfx.long.count$fxn)
head(sfx.long.count)
sfx.wide.count <- dcast(sfx.long.count, formula = full_fxn_tax ~ sampleID, value.var = "count_abun")
dim(sfx.wide.count) # 31719   102
sel.na <- which(is.na(sfx.wide.count),arr.ind = TRUE)
sfx.wide.count[sel.na] <- 0
sum(colSums(sfx.wide.count[,-1])) # 261705854
hist(colSums(sfx.wide.count[,-1]))
mean(colSums(sfx.wide.count[,-1])) # 2591147
sd(colSums(sfx.wide.count[,-1])) # 1496129

summary(colSums(sfx.wide.count[,-1]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 417218 1471726 2087559 2591147 3473391 7417074 
length(unique(sfx.long.count$subsys_L1)) # 35

fxn_sum_counts <- colSums(sfx.wide.count[,-1])

# # # #


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host-removal - get into Phyloseq object
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 31719     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 31719   101

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31719 taxa and 101 samples ]
# tax_table()   Taxonomy Table:    [ 31719 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR9276169" "SRR9276170" "SRR9276171" "SRR9276174" "SRR9276175" "SRR9276176" "SRR9276177" "SRR9276178" "SRR9276179" "SRR9276180" "SRR9276181" "SRR9276182" "SRR9276183"
# [14] "SRR9276184" "SRR9276185" "SRR9276186" "SRR9276188" "SRR9276189" "SRR9276190" "SRR9276191" "SRR9276192" "SRR9276193" "SRR9276196" "SRR9276197" "SRR9276224" "SRR9276225"
# [27] "SRR9276226" "SRR9276227" "SRR9276228" "SRR9276229" "SRR9276230" "SRR9276231" "SRR9276232" "SRR9276233" "SRR9276243" "SRR9276244" "SRR9276245" "SRR9276246" "SRR9276247"
# [40] "SRR9276248" "SRR9276249" "SRR9276251" "SRR9276252" "SRR9276253" "SRR9276254" "SRR9276255" "SRR9276256" "SRR9276257" "SRR9276258" "SRR9276259" "SRR9276260" "SRR9276269"
# [53] "SRR9276270" "SRR9276271" "SRR9276272" "SRR9276275" "SRR9276277" "SRR9276278" "SRR9276279" "SRR9276280" "SRR9276281" "SRR9276294" "SRR9276303" "SRR9276304" "SRR9276305"
# [66] "SRR9276306" "SRR9276317" "SRR9276318" "SRR9276319" "SRR9276320" "SRR9276321" "SRR9276322" "SRR9276323" "SRR9276324" "SRR9276325" "SRR9276326" "SRR9276335" "SRR9276336"
# [79] "SRR9276337" "SRR9276338" "SRR9276339" "SRR9276340" "SRR9276341" "SRR9276342" "SRR9276343" "SRR9276344" "SRR9276345" "SRR9276346" "SRR9276347" "SRR9276356" "SRR9276359"
# [92] "SRR9276360" "SRR9276361" "SRR9276362" "SRR9276363" "SRR9276364" "SRR9276365" "SRR9276366" "SRR9276367" "SRR9276368" "SRR9276369"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

head(row.names(sradat.select))

samp <- sradat.select

dim(samp) # 102  34

head(row.names(samp)) # 

row.names(samp) <- samp$Run

identical(row.names(samp), sample_names(phy)) # FALSE
length(row.names(samp)) # 102
length(sample_names(phy)) # 101
sel <- which(row.names(samp) %in% sample_names(phy)) # 101

samp2 <- samp[sample_names(phy), ]

identical(row.names(samp2), names(fxn_sum_counts)) # TRUE

samp2$fxn_sum_counts <- fxn_sum_counts


SAMP <- sample_data(samp2)



### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31719 taxa and 101 samples ]
# sample_data() Sample Data:       [ 101 samples by 35 sample variables ]
# tax_table()   Taxonomy Table:    [ 31719 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"                                        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_isomerase"                                                  
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)"
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                      
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                      
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-oxalomesaconate_tautomerase_(EC_5.3.2.8)"     


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = phy, file = "phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS")


head(phy@sam_data)

# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 101
NonZeroFxns

mean(NonZeroFxns) # 12203.89
sd(NonZeroFxns) # 2757.84


table(phy@sam_data$Group)
# D_first_sampling D_taking_environment_samples             F_first_sampling F_taking_environment_samples             G_first_sampling 
# 30                            4                           29                            4                           30 
# G_taking_environment_samples 
# 4 

# mouse fecal samples
sel <- which(phy@sam_data$Group %in% c("D_first_sampling", "G_first_sampling", "F_first_sampling" ))
mean(NonZeroFxns[sel]) # 11541.13
sd(NonZeroFxns[sel]) # 1847.001

# soil samples
mean(NonZeroFxns[-sel]) # 17119.33
sd(NonZeroFxns[-sel]) # 3455.328


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host removal - COPY of R code to run CPP steps on HPC
#    1) build reaction search - get reactions and compounds
#    2) get cpd rel abun per sample
#    3) collate compounds for each sample
#-------------------------

# # # # # # # # # # # # #
# #
# # R script for cpp3d
# # - build reaction search in parallel - get_reactions & compounds
# # - get cpd rel abun per sample
# # - collate_compounds
# #
# # For study - Liu_2020_mice_soil_PRJNA542998
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# message("\n# establish folders and input files")
# 
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis'")
# workdir <- "/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis/working"
# 
# message("\nthis_study <- '-Liu-2020-mice-pawsey'")
# this_study <- "-Liu-2020-mice-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS")
# 
# 
# subsys.lut <- readRDS("subsys.lut.RDS")
# rxns.lut <- readRDS("rxns.lut.RDS")
# rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
# compounds.lut <- readRDS("compounds.lut.RDS")
# 
# 
# 
# message("\n### 1) build reaction search in parallel - get_reactions & compounds")
# message("\n# # # # # # # # # #")
# message("\ndf.tax <- as.data.frame(phy@tax_table)")
# df.tax <- as.data.frame(phy@tax_table)
# message("\nhead(row.names(df.tax))")
# head(row.names(df.tax))
# message("\ndim(df.tax)")
# dim(df.tax)
# 
# 
# get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
#   
#   rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#   rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#   
#   sub1 <- df.tax$subsys_L1[i]
#   sub2 <- df.tax$subsys_L2[i]
#   sub3 <- df.tax$subsys_L3[i]
#   
#   fxn.temp <- df.tax$fxn[i]
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   
#   # store results corresponding to each Superfocus row
#   fxn.list <- list()
#   fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
#   
#   # check for multiple functions/reactions?
#   flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
#   flag2 <- grepl(pattern = "_@_", x = fxn.temp)
#   if (!any(flag1,flag2)==TRUE) {
#     # no multiples
#     fxns <- fxn.temp
#   } else if (flag1==TRUE) {
#     fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
#   } else {
#     fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
#   }
#   # remove underscores
#   ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
#   
#   # process each fxn & store attributes
#   df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
#   
#   # Identify '/' separators with no '_'  ??
#   
#   for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
#     #f<-1
#     #f<-2
#     f.in <- fxns[f]
#     
#     # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
#     full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
#     full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
#     
#     ## data cleaning
#     
#     # trim off '_#' and '_##' tags
#     trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
#     if (!is.na(trim_nchar) & length(trim_nchar)==1) {
#       f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
#     }
#     
#     # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
#     f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
#     
#     #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#     #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#     
#     full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
#     full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
#     
#     sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
#     
#     ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
#     if (length(sel.rx)>=1) {
#       df.fxns$matching_method[f] <- "Exact hierachy match"
#       df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#       
#     } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
#       # search by EC id if present
#       
#       f.in <- fxns[f] # this goes back to string with brackets for EC
#       ## LOOK FOR MULTIPLE ECs ?
#       
#       how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
#       
#       ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
#       #class(ECs)
#       ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
#       ECs.collapse <- paste0(ECs, collapse = "|")
#       
#       sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
#       
#       if (length(how_many_ECs)==0 | length(ECs)==0) {
#         # there was a glitch, database typo, or some error in identifying the EC number
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       } else if (length(sel.rx)>=1) {
#         # combined EC hits identified
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
#         # treat EC hits individually
#         sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
#         # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
#         sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else {
#         # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       # END EC matching
#       
#       
#     } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
#       ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
#       # 1. 'reactions' table by name: rxns.lut$name
#       # 2. 'reactions' table by aliases: rxns.lut$aliases
#       # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
#       # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
#       
#       if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
#         # 1a - exact match - rxns.lut$name
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
#         #rxns.lut$name[sel.rx]
#         df.fxns$matching_method[f] <- "Matched Reactions name"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
#         # 2a - exact match - rxns.lut$aliases
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
#         #rxns.lut$aliases[sel.rx]
#         #rxns.lut$name[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Reactions aliases"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
#         # 3a - exact match - subsys.lut$Role
#         sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
#         #subsys.lut$Role[sel.rx]
#         #subsys.lut$Reaction[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Subsytem role"
#         df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
#         # 4a - exact match - rxn_pathways.lut$External_rxn_name
#         sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
#         
#         df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
#         df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
#         
#         
#       } else {
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       
#       ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
#       
#       
#     } # END function - reaction search
#     
#     #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
#     
#     #print(paste0("completed fxn ", f))
#     
#     
#     ## now investigate these reactions ...
#     # Reactions lookup table: 
#     # - "equation": Definition of reaction expressed using compound IDs and after protonation
#     # Compounds lookup table:
#     # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
#     #df.fxns
#     
#     
#     #if (df.fxns$matching_method == "No match found") {
#     if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
#       
#       df.Rxns <- NA
#       df.Compounds <- NA
#       
#     } else { # reaction(s) were identified
#       
#       # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
#       f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
#       
#       df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
#                             rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
#       
#       for (r in 1:dim(df.Rxns)[1]) {
#         #r<-1
#         #this_rxn <- "rxn00004"
#         this_rxn <- df.Rxns$rxn_id[r]
#         sel <- which(rxns.lut$id == this_rxn)
#         ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
#         ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
#         ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
#         
#         # extract compound info
#         
#         #df.Rxns$rxn_eqn[r]
#         #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
#         #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
#         
#         ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
#         # 5 23 43 61
#         # 6 25 46 65 83
#         
#         ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
#         # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
#         
#         if (length(compds)>=1) {
#           
#           df.Rxns$compds[r] <- paste0(compds, collapse = ";")
#           
#           ## get compound coefficients?
#           start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
#           end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
#           ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
#           
#           df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
#           
#           # get formulas of compounds
#           
#           formx <-filter(compounds.lut, id %in% compds )
#           row.names(formx) <- formx$id
#           ( formx.char <- formx[compds, ]$formula )
#           # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
#           # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
#           # "C7H7O7" "H2O"    "C7H5O6"
#           df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
#           
#           ( compd.names <- formx[compds, ]$name )
#           # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
#           
#           temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
#                                           cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
#                                           #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
#           )
#           
#         } else {
#           # No specified reaction equation or chemical formula info
#           df.Rxns$compds[r] <- NA
#           df.Rxns$compd_coef[r] <- NA
#           df.Rxns$chem_formx[r] <- NA
#           
#           temp.df.Compounds <- NA
#           
#         }
#         
#         if (r==1) { df.Compounds <- temp.df.Compounds }
#         
#         if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
#         
#         # clean up - if there are additional reactions?
#         temp.df.Compounds <- NA
#         
#       } # END loop for r - rxn_id's per f/f.in
#       
#     } # END else loop when reactions identified
#     
#     # store results corresponding to each sub-reaction of each Superfocus row
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
#     
#     
#   } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
#   
#   saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
#   
# } # END function to be run in parallel for each superfocus row
# 
# 
# # # # # # # # # # # # # # # # # # #
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
#   get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# message("\n(num_results_files <- dim(df.tax)[1])")
# (num_results_files <- dim(df.tax)[1])
# 
# # assemble all compound data outputs
# # start with blank row
# 
# df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )
# 
# for (i in 1:num_results_files) {
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
#   
#   f_no <- length( temp[[1]][["compounds"]] )
#   
#   for (f in 1:f_no) {
#     #f<-2
#     # only add non-NA results
#     if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
#       
#       df.temp <- temp[[1]][["compounds"]][[f]]
#       ok <- complete.cases(df.temp)
#       df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
#       df.out <- rbind(df.out,df.temp)
#     }
#   }
#   print(paste0("added df ",i," of ",num_results_files ))
#   
# }
# 
# 
# message("\nstr(df.out)")
# str(df.out)
# 
# 
# saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))
# 
# # remove NA first row
# message("\nhead(df.out)")
# head(df.out)
# 
# df.out <- df.out[-1, ]
# 
# message("\ndim(df.out)")
# dim(df.out)
# 
# 
# message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")
# 
# df.out$cpd_molar_prop_norm <- NA
# 
# message("\nlength(unique(df.out$superfocus_fxn))")
# length(unique(df.out$superfocus_fxn))
# 
# message("\nphy")
# phy
# 
# message("\n% of functions represented - with compound information")
# 100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))
# 
# 
# fxns_found <- unique(df.out$superfocus_fxn)
# 
# for (k in 1:length(fxns_found)) {
#   #k<-1
#   this_fxn <- fxns_found[k]
#   sel <- which(df.out$superfocus_fxn == this_fxn)
#   
#   sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
#   # calculate 
#   
#   df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
#   
#   print(paste0("completed ",k))
#   
# }
# 
# message("\nsum(df.out$cpd_molar_prop_norm)")
# sum(df.out$cpd_molar_prop_norm)
# 
# message("\nsample_sums(phy)")
# sample_sums(phy)
# 
# message("\ngetwd()")
# getwd()
# 
# saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))
# 
# 
# 
# message("\n### 2) get cpd rel abun per sample")
# message("\n# # # # # # # # # #")
# 
# 
# df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
# message("\ndim(df.OTU)")
# dim(df.OTU)
# 
# 
# get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
#   
#   this_samp <- sample_names(phy_in)[i]
#   df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
#   
#   dat.cpd$sample <- this_samp
#   
#   dat.cpd$cpd_rel_abun_norm <- NA
#   
#   fxns_all <- row.names(df.OTU)
#   
#   for (k in 1:length(fxns_all)) {
#     #k<-1
#     this_fxn <- fxns_all[k]
#     sel <- which(dat.cpd$superfocus_fxn == this_fxn)
#     
#     if (length(sel)>=1) {
#       dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
#       
#     }
#   } # END rel abun values for all relevant functions added
#   
#   saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   
# } # END
# 
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
# head(dat)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   dat <- rbind(dat, temp)
#   print(paste0("completed ",i))
# }
# 
# 
# saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )
# 
# rm(temp)
# 
# message("\nstr(dat)")
# str(dat)
# 
# message("\nsum(dat$cpd_rel_abun_norm)")
# sum(dat$cpd_rel_abun_norm)
# 
# message("\naverage functional relative abundance per sample")
# message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# 
# message("\nnames(dat)")
# names(dat)
# 
# message("\nlength(unique(dat$cpd_id))")
# length(unique(dat$cpd_id))
# 
# 
# 
# 
# message("\n### 3) collate_compounds within each sample")
# message("\n# # # # # # # # # #")
# 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
# } # END
# 
# 
# 
# no_forks <- 4
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# message("\nstr(dat.cpd.collate)")
# str(dat.cpd.collate)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)")
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))
# 
# # END


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC
#-------------------------

# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "4.1"
# 
# $year
# [1] "2024"
# 
# $month
# [1] "06"
# 
# $day
# [1] "14"
# 
# $`svn rev`
# [1] "86737"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.4.1 (2024-06-14)"
# 
# $nickname
# [1] "Race for Your Life"
# 
# [1] ‘4.4.1’
# Loading required package: foreach
# Loading required package: iterators
# [1] ‘1.0.17’
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# [1] ‘1.1.4’
# [1] ‘1.5.2’
# [1] ‘1.46.0’
# 
# # establish folders and input files
# 
# workdir <- '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis/working'
# 
# this_study <- '-Liu-2020-mice-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS')
# 
# ### 1) build reaction search in parallel - get_reactions & compounds
# 
# # # # # # # # # # #
# 
# df.tax <- as.data.frame(phy@tax_table)
# 
# head(row.names(df.tax))
# [1] "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"
# 
# dim(df.tax)
# [1] 31719     4
# [[1]]
# NULL
# 
# ...
# 
# 
# [[31719]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 31719
# [1] "added df 1 of 31719"
# [1] "added df 2 of 31719"
# [1] "added df 3 of 31719"
# ...
# 
# 
# [1] "added df 31717 of 31719"
# [1] "added df 31718 of 31719"
# [1] "added df 31719 of 31719"
# 
# str(df.out)
# 'data.frame':	1108221 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id        : chr  NA "rxn25278" "rxn25278" "rxn25279" ...
# $ cpd_id        : chr  NA "cpd25681" "cpd02597" "cpd24620" ...
# $ cpd_name      : chr  NA "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form      : chr  NA "C7H5O6" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f
# 1           <NA> NA
# 2          fxn_2  1
# 3          fxn_2  1
# 4          fxn_3  1
# 5          fxn_3  1
# 6          fxn_3  1
# f__in
# 1                                                                          <NA>
#   2                                                   2-methylaconitate isomerase
# 3                                                   2-methylaconitate isomerase
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 6 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn25278 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 3 rxn25278 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 4 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 5 rxn25279 cpd00001                                          H2O      H2O
# 6 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 1108220       8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 16756
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31719 taxa and 101 samples ]
# sample_data() Sample Data:       [ 101 samples by 35 sample variables ]
# tax_table()   Taxonomy Table:    [ 31719 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 52.82638
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 16755"
# [1] "completed 16756"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 16756
# 
# sample_sums(phy)
# SRR9276169 SRR9276170 SRR9276171 SRR9276174 SRR9276175 SRR9276176 SRR9276177 
# 100        100        100        100        100        100        100 
# SRR9276178 SRR9276179 SRR9276180 SRR9276181 SRR9276182 SRR9276183 SRR9276184 
# 100        100        100        100        100        100        100 
# SRR9276185 SRR9276186 SRR9276188 SRR9276189 SRR9276190 SRR9276191 SRR9276192 
# 100        100        100        100        100        100        100 
# SRR9276193 SRR9276196 SRR9276197 SRR9276224 SRR9276225 SRR9276226 SRR9276227 
# 100        100        100        100        100        100        100 
# SRR9276228 SRR9276229 SRR9276230 SRR9276231 SRR9276232 SRR9276233 SRR9276243 
# 100        100        100        100        100        100        100 
# SRR9276244 SRR9276245 SRR9276246 SRR9276247 SRR9276248 SRR9276249 SRR9276251 
# 100        100        100        100        100        100        100 
# SRR9276252 SRR9276253 SRR9276254 SRR9276255 SRR9276256 SRR9276257 SRR9276258 
# 100        100        100        100        100        100        100 
# SRR9276259 SRR9276260 SRR9276269 SRR9276270 SRR9276271 SRR9276272 SRR9276275 
# 100        100        100        100        100        100        100 
# SRR9276277 SRR9276278 SRR9276279 SRR9276280 SRR9276281 SRR9276294 SRR9276303 
# 100        100        100        100        100        100        100 
# SRR9276304 SRR9276305 SRR9276306 SRR9276317 SRR9276318 SRR9276319 SRR9276320 
# 100        100        100        100        100        100        100 
# SRR9276321 SRR9276322 SRR9276323 SRR9276324 SRR9276325 SRR9276326 SRR9276335 
# 100        100        100        100        100        100        100 
# SRR9276336 SRR9276337 SRR9276338 SRR9276339 SRR9276340 SRR9276341 SRR9276342 
# 100        100        100        100        100        100        100 
# SRR9276343 SRR9276344 SRR9276345 SRR9276346 SRR9276347 SRR9276356 SRR9276359 
# 100        100        100        100        100        100        100 
# SRR9276360 SRR9276361 SRR9276362 SRR9276363 SRR9276364 SRR9276365 SRR9276366 
# 100        100        100        100        100        100        100 
# SRR9276367 SRR9276368 SRR9276369 
# 100        100        100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/Liu_2020_mice_soil_PRJNA542998/cpp_analysis"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 31719   101
# [[1]]
# NULL
# 
# ...
# 
# 
# 
# [[101]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f
# 2          fxn_2 1
# 3          fxn_2 1
# 4          fxn_3 1
# 5          fxn_3 1
# 6          fxn_3 1
# 7          fxn_4 1
# f__in
# 2                                                   2-methylaconitate isomerase
# 3                                                   2-methylaconitate isomerase
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 6 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 7                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn25278 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 3 rxn25278 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 4 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 5 rxn25279 cpd00001                                          H2O      H2O
# 6 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 7 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# cpd_molar_prop cpd_molar_prop_norm     sample cpd_rel_abun_norm
# 2              1          0.50000000 SRR9276169      0.000000e+00
# 3              1          0.50000000 SRR9276169      0.000000e+00
# 4              1          0.33333333 SRR9276169      7.421797e-06
# 5              1          0.33333333 SRR9276169      7.421797e-06
# 6              1          0.33333333 SRR9276169      7.421797e-06
# 7              1          0.05555556 SRR9276169      0.000000e+00
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# ...
# 
# 
# [1] "completed 100"
# [1] "completed 101"
# 
# str(dat)
# 'data.frame':	111930220 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...
# $ sample             : chr  "SRR9276169" "SRR9276169" "SRR9276169" "SRR9276169" ...
# $ cpd_rel_abun_norm  : num  0.00 0.00 7.42e-06 7.42e-06 7.42e-06 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 6926.587
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 68.58007
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7967
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# [[1]]
# NULL
# 
# ...
# 
# 
# [[101]]
# NULL
# 
# 
# ## assemble results
# cpd_id     sample cpd_rel_abun
# 1 cpd25681 SRR9276169 6.008801e-05
# 2 cpd02597 SRR9276169 7.466500e-02
# 3 cpd24620 SRR9276169 5.213608e-05
# 4 cpd00001 SRR9276169 4.880965e+00
# 5 cpd01501 SRR9276169 4.842088e-02
# 6 cpd00851 SRR9276169 2.354891e-02
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# ...
# 
# 
# 
# [1] "completed 100"
# [1] "completed 101"
# 
# str(dat.cpd.collate)
# 'data.frame':	804667 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9276169" "SRR9276169" "SRR9276169" "SRR9276169" ...
# $ cpd_rel_abun: num  6.01e-05 7.47e-02 5.21e-05 4.88 4.84e-02 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 6926.587
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 68.58007
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host-removal - continue CPP analysis
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Liu_2020_mice_soil_PRJNA542998/cpp_analysis/dat.cpd.collate-all-samps-cpp3d--Liu-2020-mice-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	804667 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9276169" "SRR9276169" "SRR9276169" "SRR9276169" ...
# $ cpd_rel_abun: num  6.01e-05 7.47e-02 5.21e-05 4.88 4.84e-02 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000068 0.008608 0.000932 7.891812 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -6.4427 -4.1692    -Inf -3.0306  0.8972 


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -9.4328 -6.4427 -4.1692 -4.9522 -3.0306  0.8972 


# make group variable from sample name

dat.cpd.collate$group <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31719 taxa and 101 samples ]
# sample_data() Sample Data:       [ 101 samples by 35 sample variables ]
# tax_table()   Taxonomy Table:    [ 31719 taxa by 4 taxonomic ranks ]

samp <- as(phy@sam_data, "data.frame")

unique(samp$Group)
# [1] "G_taking_environment_samples" "D_taking_environment_samples" "F_first_sampling"             "F_taking_environment_samples" "D_first_sampling"            
# [6] "G_first_sampling"   


for (i in 1:length(sample_names(phy))) {
  #for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  #this_samp <- samp$Run[i]
  sel.phy <- which(phy@sam_data$Run == this_samp)
  sel.dat <- which(dat.cpd.collate$sample == this_samp)
  
  dat.cpd.collate$group[sel.dat] <- phy@sam_data$Group[sel.phy]
  
  #dat.cpd.collate$group[sel.dat] <- as.character( samp$group_new[i] )
  
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) #
# [1] "G_taking_environment_samples" "D_taking_environment_samples" "F_first_sampling"             "F_taking_environment_samples" "D_first_sampling"            
# [6] "G_first_sampling" 

dat.cpd.collate$group_label <- factor(dat.cpd.collate$group, 
                                      levels = c("D_first_sampling", "G_first_sampling", "F_first_sampling",
                                                 "D_taking_environment_samples", "G_taking_environment_samples", "F_taking_environment_samples" ),
                                      labels = c("Desert", "Grassland", "Forest",
                                                 "Desert soil", "Grassland soil", "Forest soil" ),
                                      ordered = TRUE)

head(dat.cpd.collate)
# cpd_id     sample cpd_rel_abun log10_abun                        group    group_label
# 1 cpd25681 SRR9276169 6.008801e-05 -4.2212122 G_taking_environment_samples Grassland soil
# 2 cpd02597 SRR9276169 7.466500e-02 -1.1268829 G_taking_environment_samples Grassland soil
# 3 cpd24620 SRR9276169 5.213608e-05 -4.2828616 G_taking_environment_samples Grassland soil
# 4 cpd00001 SRR9276169 4.880965e+00  0.6885057 G_taking_environment_samples Grassland soil
# 5 cpd01501 SRR9276169 4.842088e-02 -1.3149674 G_taking_environment_samples Grassland soil
# 6 cpd00851 SRR9276169 2.354891e-02 -1.6280291 G_taking_environment_samples Grassland soil

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Liu_2020_mice_soil_PRJNA542998-Hostremoval.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Liu_2020_mice_soil_PRJNA542998-Hostremoval.rds")


str(dat.cpd.collate)
# 'data.frame':	804667 obs. of  6 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9276169" "SRR9276169" "SRR9276169" "SRR9276169" ...
# $ cpd_rel_abun: num  6.01e-05 7.47e-02 5.21e-05 4.88 4.84e-02 ...
# $ log10_abun  : num  -4.221 -1.127 -4.283 0.689 -1.315 ...
# $ group       : chr  "G_taking_environment_samples" "G_taking_environment_samples" "G_taking_environment_samples" "G_taking_environment_samples" ...
# $ group_label : Ord.factor w/ 6 levels "Desert"<"Grassland"<..: 5 5 5 5 5 5 5 5 5 5 ...


length( unique(dat.cpd.collate$cpd_id) ) # 7967
7967*101 # 804667


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id     sample cpd_rel_abun log10_abun                        group    group_label
# 1 cpd25681 SRR9276169 6.008801e-05 -4.2212122 G_taking_environment_samples Grassland soil
# 2 cpd02597 SRR9276169 7.466500e-02 -1.1268829 G_taking_environment_samples Grassland soil
# 3 cpd24620 SRR9276169 5.213608e-05 -4.2828616 G_taking_environment_samples Grassland soil
# 4 cpd00001 SRR9276169 4.880965e+00  0.6885057 G_taking_environment_samples Grassland soil
# 5 cpd01501 SRR9276169 4.842088e-02 -1.3149674 G_taking_environment_samples Grassland soil
# 6 cpd00851 SRR9276169 2.354891e-02 -1.6280291 G_taking_environment_samples Grassland soil

dim(data_in) # 804667      6

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 6481.139
sd(no_compounds) #  452.119

mean(sample_sum_relabun) # 68.58007
sd(sample_sum_relabun) # 2.952826

#length(unique(data_in$cpd_id)) # 7967
length(unique(data_in$cpd_id[ which(data_in$cpd_rel_abun > 0) ])) # 7967


## Only gut samples

data_in <- filter( dat.cpd.collate, group_label %in% c("Desert","Grassland","Forest") )

dim(data_in) # 709063      6
unique_samps <- unique(data_in$sample)
no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 6381.629
sd(no_compounds) #  355.4287

mean(sample_sum_relabun) # 68.49007
sd(sample_sum_relabun) # 3.106018

length(unique(data_in$cpd_id[ which(data_in$cpd_rel_abun > 0) ])) # 7782


## Only soil samples

data_in <- filter( dat.cpd.collate, group_label %in% c("Desert soil","Grassland soil","Forest soil") )

dim(data_in) # 95604     6
unique_samps <- unique(data_in$sample)
no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 7219.167
sd(no_compounds) #  416.3309

mean(sample_sum_relabun) # 69.24758
sd(sample_sum_relabun) # 1.2389

length(unique(data_in$cpd_id[ which(data_in$cpd_rel_abun > 0) ])) # 7913


#-------------------------

#### Liu_2020_mice_soil_PRJNA542998 - Mouse environment exposures - w/ Host-removal - continue CPP analysis
#    CPP - get into phyloseq object
#    beta diversity
#    alpha diversity
#    response in selected compounds: glucose, cellulose, CO2, O2, AEC, ATP/ADP
#    heatmap of scaled CPP
#-------------------------

phy.fxn <- readRDS('phy-phyloseq-fxn-Liu_2020_mice_soil_PRJNA542998-Host-removal.RDS')
phy.fxn
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 31719 taxa and 101 samples ]
# sample_data() Sample Data:       [ 101 samples by 35 sample variables ]
# tax_table()   Taxonomy Table:    [ 31719 taxa by 4 taxonomic ranks ]

table(phy.fxn@sam_data$Group)
# D_first_sampling D_taking_environment_samples             F_first_sampling F_taking_environment_samples 
# 30                            4                           29                            4 
# G_first_sampling G_taking_environment_samples 
# 30                            4 

str(phy.fxn@sam_data)
# 'data.frame':	101 obs. of  35 variables:
#   Formal class 'sample_data' [package "phyloseq"] with 4 slots
# ..@ .Data    :List of 35
# .. ..$ : chr  "SRR9276169" "SRR9276170" "SRR9276171" "SRR9276174" ...
# .. ..$ : chr  "WGS" "WGS" "WGS" "WGS" ...
# .. ..$ : num  286 264 270 252 252 252 252 252 252 252 ...
# .. ..$ : num  1.24e+10 1.48e+10 1.27e+10 6.83e+09 6.86e+09 ...
# .. ..$ : chr  "PRJNA542998" "PRJNA542998" "PRJNA542998" "PRJNA542998" ...
# .. ..$ : chr  "SAMN11664238" "SAMN11664231" "SAMN11664232" "SAMN11664340" ...
# .. ..$ : chr  "Metagenome or environmental" "Metagenome or environmental" "Metagenome or environmental" "Metagenome or environmental" ...
# .. ..$ : num  4.94e+09 6.23e+09 5.24e+09 2.68e+09 2.69e+09 ...
# .. ..$ : chr  "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" "INNER MONGOLIA AGRICULTURAL UNIVERSITY" ...
# .. ..$ : POSIXct, format: "2014-12-01" "2014-12-01" "2014-12-01" "2014-12-01" ...
# .. ..$ : chr  "public" "public" "public" "public" ...
# .. ..$ : chr  "sra,run.zq,fastq" "run.zq,fastq,sra" "run.zq,fastq,sra" "fastq,sra,run.zq" ...
# .. ..$ : chr  "ncbi,s3,gs" "s3,gs,ncbi" "ncbi,gs,s3" "gs,s3,ncbi" ...
# .. ..$ : chr  "ncbi.public,s3.us-east-1,gs.us-east1" "gs.us-east1,s3.us-east-1,ncbi.public" "s3.us-east-1,ncbi.public,gs.us-east1" "s3.us-east-1,gs.us-east1,ncbi.public" ...
# .. ..$ : chr  "SRX6046150" "SRX6046149" "SRX6046148" "SRX6046145" ...
# .. ..$ : chr  "China" "China" "China" "China" ...
# .. ..$ : chr  "Asia" "Asia" "Asia" "Asia" ...
# .. ..$ : chr  "China" "China" "China" "China" ...
# .. ..$ : chr  "mice" "mice" "mice" "mice" ...
# .. ..$ : chr  "Illumina HiSeq 2500" "Illumina HiSeq 2500" "Illumina HiSeq 2500" "Illumina HiSeq 2500" ...
# .. ..$ : chr  "mice fecals" "mice fecals" "mice fecals" "mice fecals" ...
# .. ..$ : chr  "29.35 N 106.33 E" "29.35 N 106.33 E" "29.35 N 106.33 E" "29.35 N 106.33 E" ...
# .. ..$ : chr  "D13247" "D13240" "D13241" "D13349" ...
# .. ..$ : chr  "PAIRED" "PAIRED" "PAIRED" "PAIRED" ...
# .. ..$ : chr  "DNase" "DNase" "DNase" "DNase" ...
# .. ..$ : chr  "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" "METAGENOMIC" ...
# .. ..$ : chr  "feces metagenome" "feces metagenome" "feces metagenome" "feces metagenome" ...
# .. ..$ : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
# .. ..$ : chr  "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" "2020-06-14T00:00:00Z" ...
# .. ..$ : chr  "2019-06-11T22:08:00Z" "2019-06-11T22:12:00Z" "2019-06-11T22:09:00Z" "2019-06-11T22:00:00Z" ...
# .. ..$ : num  1 1 1 1 1 1 1 1 1 1 ...
# .. ..$ : chr  "D13247" "D13240" "D13241" "D13349" ...
# .. ..$ : chr  "SRP201145" "SRP201145" "SRP201145" "SRP201145" ...
# .. ..$ : chr  "G_taking_environment_samples" "D_taking_environment_samples" "D_taking_environment_samples" "F_first_sampling" ...
# .. ..$ : num  4191857 5618853 4887846 1392392 1471726 ...
# ..@ names    : chr  "Run" "Assay Type" "AvgSpotLen" "Bases" ...
# ..@ row.names: chr  "SRR9276169" "SRR9276170" "SRR9276171" "SRR9276174" ...
# ..@ .S3Class : chr "data.frame"


phy.fxn <- prune_samples( samples = phy.fxn@sam_data$Group %in% c("D_first_sampling", "G_first_sampling", "F_first_sampling"), x = phy.fxn )

phy.fxn@sam_data$group_label <- factor( phy.fxn@sam_data$Group, levels = c("D_first_sampling", "G_first_sampling", "F_first_sampling"),
                                        labels = c("Desert", "Grassland", "Forest"), ordered = TRUE)
samp <- as( phy.fxn@sam_data, "data.frame")
temp <- samp


dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Liu_2020_mice_soil_PRJNA542998-Hostremoval.rds")

str(dat.cpd.collate)
# 'data.frame':	804667 obs. of  6 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9276169" "SRR9276169" "SRR9276169" "SRR9276169" ...
# $ cpd_rel_abun: num  6.01e-05 7.47e-02 5.21e-05 4.88 4.84e-02 ...
# $ log10_abun  : num  -4.221 -1.127 -4.283 0.689 -1.315 ...
# $ group       : chr  "G_taking_environment_samples" "G_taking_environment_samples" "G_taking_environment_samples" "G_taking_environment_samples" ...
# $ group_label : Ord.factor w/ 6 levels "Desert"<"Grassland"<..: 5 5 5 5 5 5 5 5 5 5 ...

#data_in <- dat.cpd.collate


data_in <- filter( dat.cpd.collate, group_label %in% c("Desert","Grassland","Forest") )

str(data_in)
# 'data.frame':	709063 obs. of  6 variables:

length( unique(data_in$cpd_id) ) # 7967
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 7782
length( unique(data_in$sample) ) # 89



### get data into phyloseq object ...

head(data_in)
#     cpd_id     sample cpd_rel_abun log10_abun            group group_label
# 1 cpd25681 SRR9276174 2.633358e-05 -4.5794901 F_first_sampling      Forest
# 2 cpd02597 SRR9276174 2.214690e-04 -3.6546871 F_first_sampling      Forest
# 3 cpd24620 SRR9276174 2.633358e-05 -4.5794901 F_first_sampling      Forest
# 4 cpd00001 SRR9276174 4.882249e+00  0.6886199 F_first_sampling      Forest
# 5 cpd01501 SRR9276174 2.074767e-04 -3.6830307 F_first_sampling      Forest
# 6 cpd00851 SRR9276174 1.344470e-03 -2.8714489 F_first_sampling      Forest


df.wide <- dcast(data_in, formula = sample + group_label ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#       sample group_label cpd00001 cpd00002  cpd00003  cpd00004  cpd00005  cpd00006   cpd00007 cpd00008
# 1 SRR9276174      Forest 4.882249 3.263526 0.6423462 0.6192881 0.5917658 0.5919520 0.04660319 2.046516
# 2 SRR9276175      Forest 4.891980 3.222542 0.5964531 0.5734767 0.5423694 0.5426748 0.05500806 1.991837
# 3 SRR9276176      Forest 4.904033 3.264556 0.6522388 0.6289964 0.5935464 0.5936672 0.04470705 2.038716
# 4 SRR9276177      Forest 4.676776 3.239831 0.6381518 0.6130299 0.5670377 0.5673722 0.05560286 1.999844
# 5 SRR9276178      Forest 4.924885 3.251344 0.6420437 0.6202340 0.5952900 0.5952373 0.04396980 2.022077

# save group variable
samp <- df.wide[ ,1:2]
#samp <- df.wide[ ,1:3]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group' column
#df.wide <- t(df.wide[ ,-c(2,3)]) # minus 'Organism' and 'abbrev' columns

head(df.wide)

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array"
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  7967   89
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix.
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)


# # convert Taxonomy table to matrix

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7967 taxa and 89 samples ]
# tax_table()   Taxonomy Table:    [ 7967 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# [1] "SRR9276174" "SRR9276175" "SRR9276176" "SRR9276177" "SRR9276178" "SRR9276179" "SRR9276180" "SRR9276181"
# [9] "SRR9276182" "SRR9276183" "SRR9276184" "SRR9276185" "SRR9276186" "SRR9276188" "SRR9276189" "SRR9276190"
# [17] "SRR9276191" "SRR9276192" "SRR9276193" "SRR9276224" "SRR9276225" "SRR9276226" "SRR9276227" "SRR9276228"
# [25] "SRR9276229" "SRR9276230" "SRR9276231" "SRR9276233" "SRR9276243" "SRR9276244" "SRR9276245" "SRR9276246"
# [33] "SRR9276247" "SRR9276248" "SRR9276249" "SRR9276251" "SRR9276252" "SRR9276253" "SRR9276254" "SRR9276255"
# [41] "SRR9276256" "SRR9276257" "SRR9276258" "SRR9276259" "SRR9276260" "SRR9276269" "SRR9276270" "SRR9276271"
# [49] "SRR9276272" "SRR9276275" "SRR9276277" "SRR9276278" "SRR9276279" "SRR9276280" "SRR9276294" "SRR9276303"
# [57] "SRR9276304" "SRR9276305" "SRR9276306" "SRR9276317" "SRR9276318" "SRR9276319" "SRR9276320" "SRR9276321"
# [65] "SRR9276322" "SRR9276323" "SRR9276324" "SRR9276325" "SRR9276326" "SRR9276337" "SRR9276339" "SRR9276340"
# [73] "SRR9276341" "SRR9276344" "SRR9276345" "SRR9276346" "SRR9276347" "SRR9276356" "SRR9276359" "SRR9276360"
# [81] "SRR9276361" "SRR9276362" "SRR9276363" "SRR9276364" "SRR9276365" "SRR9276366" "SRR9276367" "SRR9276368"
# [89] "SRR9276369"

#identical(sample_names(phy.cpp), samp$sample) # TRUE
#identical(sample_names(phy.cpp), sradat.select2$Run) # TRUE
identical(sample_names(phy.cpp), temp$Run) # TRUE

#row.names(sradat.select2) <- sradat.select2$Run
identical( row.names(temp), sample_names(phy.cpp) ) # TRUE

#samp <- sradat.select2
samp <- temp

# row.names need to match sample_names() from phyloseq object
#row.names(samp) <- samp$sample
#identical(row.names(samp), samp$sample) # TRUE



### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp

# check for 'taxa' (compounds) with zero data - because soil samples were excluded
min(taxa_sums(phy.cpp)) # 0
# prune taxa that have zero sequence reads
phy.cpp <- prune_taxa(taxa = taxa_sums(phy.cpp) > 0, x = phy.cpp)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7782 taxa and 89 samples ]
# sample_data() Sample Data:       [ 89 samples by 36 sample variables ]
# tax_table()   Taxonomy Table:    [ 7782 taxa by 1 taxonomic ranks ]

phy.cpp@sam_data

min(taxa_sums(phy.cpp)) #  9.833847e-10

saveRDS(object = phy.cpp, file = "phy.cpp-cleaned-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.RDS")
phy.cpp <- readRDS("phy.cpp-cleaned-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.RDS")


phy_in <- phy.cpp

sum(sample_sums(phy_in)) # 6095.616
sample_sums(phy_in)
# SRR9276174 SRR9276175 SRR9276176 SRR9276177 SRR9276178 SRR9276179 SRR9276180 SRR9276181 SRR9276182 SRR9276183 
# 71.74148   70.73691   71.81279   70.82755   71.90445   71.57646   71.71403   71.38909   71.81819   72.23729 
# SRR9276184 SRR9276185 SRR9276186 SRR9276188 SRR9276189 SRR9276190 SRR9276191 SRR9276192 SRR9276193 SRR9276224 
# 70.30826   71.01337   71.08037   71.01778   72.57866   69.22789   72.13741   68.93393   71.28015   66.05282 
# SRR9276225 SRR9276226 SRR9276227 SRR9276228 SRR9276229 SRR9276230 SRR9276231 SRR9276233 SRR9276243 SRR9276244 
# 64.18318   64.85526   68.26635   67.15102   63.84957   61.04372   66.67904   71.47422   70.15545   70.82121 
# SRR9276245 SRR9276246 SRR9276247 SRR9276248 SRR9276249 SRR9276251 SRR9276252 SRR9276253 SRR9276254 SRR9276255 
# 64.18512   70.14218   70.77799   63.25889   66.39699   71.04983   68.13160   70.25744   71.00407   70.65403 
# SRR9276256 SRR9276257 SRR9276258 SRR9276259 SRR9276260 SRR9276269 SRR9276270 SRR9276271 SRR9276272 SRR9276275 
# 70.63076   71.22225   71.58401   70.89272   71.39949   70.46127   71.09305   70.98812   70.26220   66.34135 
# SRR9276277 SRR9276278 SRR9276279 SRR9276280 SRR9276294 SRR9276303 SRR9276304 SRR9276305 SRR9276306 SRR9276317 
# 66.57356   66.14565   64.82085   68.26827   56.33603   65.62894   65.77766   66.10524   66.11609   71.54200 
# SRR9276318 SRR9276319 SRR9276320 SRR9276321 SRR9276322 SRR9276323 SRR9276324 SRR9276325 SRR9276326 SRR9276337 
# 71.56440   70.97877   69.55837   70.92634   63.91269   69.57655   68.13116   71.38427   71.45713   65.58748 
# SRR9276339 SRR9276340 SRR9276341 SRR9276344 SRR9276345 SRR9276346 SRR9276347 SRR9276356 SRR9276359 SRR9276360 
# 66.51634   66.37154   65.55667   65.35271   67.86952   63.71659   69.93982   71.00578   64.64150   65.70813 
# SRR9276361 SRR9276362 SRR9276363 SRR9276364 SRR9276365 SRR9276366 SRR9276367 SRR9276368 SRR9276369 
# 67.17846   66.07311   66.65982   66.71817   66.05269   66.36931   63.66604   66.61459   70.61257

summary( sample_sums(phy_in) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 56.34   66.12   69.58   68.49   71.02   72.58 

sd( sample_sums(phy_in) )
# 3.106018

max(taxa_sums(phy_in)) # 601.1415


# don't rarefy - already in form of relative abundance %


## ordination plot
## PCoA + Bray-Curtis


## Already normalised to relative abundance

# # rarefy #1
# seed <- 1234
# r1 <- rarefy_even_depth(phy_in, sample.size = min(sample_sums(phy_in)),
#                         rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# min(taxa_sums(r1)) # 1
# sample_sums(r1) # all 3073
# ntaxa(r1) # 1014

r1 <- phy_in


### ORDINATION PLOT # # # # # # # # # # # # # # # 
### PCoA + Bray-Curtis

set.seed(1234)
ord <- ordinate(r1, "PCoA", "bray")

ord

str(r1@sam_data)

names(r1@sam_data)
# [1] "Run"                            "Assay.Type"                     "AvgSpotLen"                    
# [4] "Bases"                          "BioProject"                     "BioSample"                     
# [7] "BioSampleModel"                 "Bytes"                          "Center.Name"                   
# [10] "Collection_Date"                "Consent"                        "DATASTORE.filetype"            
# [13] "DATASTORE.provider"             "DATASTORE.region"               "Experiment"                    
# [16] "geo_loc_name_country"           "geo_loc_name_country_continent" "geo_loc_name"                  
# [19] "HOST"                           "Instrument"                     "isolation_source"              
# [22] "lat_lon"                        "Library.Name"                   "LibraryLayout"                 
# [25] "LibrarySelection"               "LibrarySource"                  "Organism"                      
# [28] "Platform"                       "ReleaseDate"                    "create_date"                   
# [31] "version"                        "Sample.Name"                    "SRA.Study"                     
# [34] "Group"                          "fxn_sum_counts"                 "group_label"  




saveRDS(r1, file = "r1-cpp-phyloseq-object-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.RDS")


p <- plot_ordination(r1, ord, type="samples", color="group_label")
#p <- plot_ordination(r1, ord, type="samples", color="Treatment_no_description", shape = "Treatment_no_description")
p

str(p$data)

# x_lab <- p$labels$x
# y_lab <- p$labels$y

x_lab <- gsub(pattern = "Axis.1", replacement = "PCo1" , x = p$labels$x)
y_lab <- gsub(pattern = "Axis.2", replacement = "PCo2" , x =  p$labels$y)


names(p$data)
# [1] "Axis.1"                         "Axis.2"                         "Run"                           
# [4] "Assay.Type"                     "AvgSpotLen"                     "Bases"                         
# [7] "BioProject"                     "BioSample"                      "BioSampleModel"                
# [10] "Bytes"                          "Center.Name"                    "Collection_Date"               
# [13] "Consent"                        "DATASTORE.filetype"             "DATASTORE.provider"            
# [16] "DATASTORE.region"               "Experiment"                     "geo_loc_name_country"          
# [19] "geo_loc_name_country_continent" "geo_loc_name"                   "HOST"                          
# [22] "Instrument"                     "isolation_source"               "lat_lon"                       
# [25] "Library.Name"                   "LibraryLayout"                  "LibrarySelection"              
# [28] "LibrarySource"                  "Organism"                       "Platform"                      
# [31] "ReleaseDate"                    "create_date"                    "version"                       
# [34] "Sample.Name"                    "SRA.Study"                      "Group"                         
# [37] "fxn_sum_counts"                 "group_label" 


cols.ecosystem_type <- c("Desert" = "gold",
                         "Grassland" = "springgreen4",
                         "Forest" = "purple4")

label = "PERMANOVA: ~Group\nR^2 = 0.57, P = 0.001\nBeta-dispersion: P = 0.36"

pp <- ggplot(data=p$data, aes(x=Axis.1, y=Axis.2)) + # , colour=Sample_type__row_type x=NMDS1, y=NMDS2
  theme_bw() + 
  
  #geom_point(aes(colour=abbrev), size = 2) + # , alpha = 0.6
  #scale_color_manual(values=cols.bact, name = "Bacteria")+
  geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  #geom_text_repel(aes(label = abbrev), size = 3)+
  geom_text_npc(npcx = "left", npcy = "top", label = label, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE, 
  
  xlab(x_lab) + ylab(y_lab)+
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.title = element_text(size = rel(1.1)),
    
    legend.background = element_rect(fill = "transparent"),
    #legend.margin = margin(t = 0,r = 0,b = 0,l = 1,unit = "pt"),
    legend.margin = margin(t = 0,r = 0,b = 0,l = -5,unit = "pt"),
    legend.key.size = unit(0.8,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8))
  )

pp

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Beta-diversity-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 10, units = "cm", res = 600, compression = "lzw",type="cairo")



## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(r1, method = "bray")
sampledf <- data.frame(sample_data(r1))
str(sampledf)


names(r1@sam_data)
# [1] "Run"                            "Assay.Type"                     "AvgSpotLen"                    
# [4] "Bases"                          "BioProject"                     "BioSample"                     
# [7] "BioSampleModel"                 "Bytes"                          "Center.Name"                   
# [10] "Collection_Date"                "Consent"                        "DATASTORE.filetype"            
# [13] "DATASTORE.provider"             "DATASTORE.region"               "Experiment"                    
# [16] "geo_loc_name_country"           "geo_loc_name_country_continent" "geo_loc_name"                  
# [19] "HOST"                           "Instrument"                     "isolation_source"              
# [22] "lat_lon"                        "Library.Name"                   "LibraryLayout"                 
# [25] "LibrarySelection"               "LibrarySource"                  "Organism"                      
# [28] "Platform"                       "ReleaseDate"                    "create_date"                   
# [31] "version"                        "Sample.Name"                    "SRA.Study"                     
# [34] "Group"                          "fxn_sum_counts"                 "group_label" 


# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df SumOfSqs      R2      F Pr(>F)    
# group_label  2  0.30933 0.56815 56.571  0.001 ***
# Residual    86  0.23513 0.43185                  
# Total       88  0.54446 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     2 0.002655 0.0013273 1.0932    999  0.362
# Residuals 86 0.104419 0.0012142 




### Alpha diversity

# Shannon Diversity Index
a_div <- plot_richness(r1, measures=c("Shannon")) #, "Simpson")) # Observed = Richness, but requires count data
a_div
# Shannon index emphasises richness, while Simpson index emphasises evenness

names(a_div$data)
# [1] "Run"                            "Assay.Type"                     "AvgSpotLen"                    
# [4] "Bases"                          "BioProject"                     "BioSample"                     
# [7] "BioSampleModel"                 "Bytes"                          "Center.Name"                   
# [10] "Collection_Date"                "Consent"                        "DATASTORE.filetype"            
# [13] "DATASTORE.provider"             "DATASTORE.region"               "Experiment"                    
# [16] "geo_loc_name_country"           "geo_loc_name_country_continent" "geo_loc_name"                  
# [19] "HOST"                           "Instrument"                     "isolation_source"              
# [22] "lat_lon"                        "Library.Name"                   "LibraryLayout"                 
# [25] "LibrarySelection"               "LibrarySource"                  "Organism"                      
# [28] "Platform"                       "ReleaseDate"                    "create_date"                   
# [31] "version"                        "Sample.Name"                    "SRA.Study"                     
# [34] "Group"                          "fxn_sum_counts"                 "group_label"                   
# [37] "samples"                        "variable"                       "value"                         
# [40] "se"  


head(a_div$data)

# include Kruskal-Wallis test results - per below: 
#ktresult <- paste0("Kruskal-Wallis\nP = 7.6 x ",expression(10^-10))
#ktresult <- "Kruskal-Wallis\nP == 7.6 * 10^-10"
ktresult <- paste0("~Kruskal-Wallis:","~P == 7.6","~x ","~10^{-10}") # "\n",
#val <- "Kruskal-Wallis\nP = 7.6 x"
#label_expr <- bquote(.(val)~10^-10)

set.seed(123)
p <- ggplot(data=a_div$data, aes(x=group_label, y=value)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  #annotate("text_npc", npcx = "right", npcy = "top", label = deparse(label_expr), parse = TRUE, size = 3 , lineheight = 0.85 )+
  annotate("text", x = c(1.2, 2.2, 3.2), y = c(5.5, 5.45, 5.45), label = c("a","b","b"), size = 3.75)+
  
  labs(x = NULL, y = "Shannon diversity") +
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Alpha-diversity-Shannon-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


names(a_div$data)
unique(a_div$data$variable) # "Shannon"

ktdat <- filter(a_div$data[ ,c("variable","value","group_label")], variable == "Shannon")

# Kruskal-Wallis test
kt <- kruskal.test( value ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  value by group_label
# Kruskal-Wallis chi-squared = 41.996, df = 2, p-value = 7.597e-10

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector

unique( ktdat$group_label )
# [1] Forest    Desert    Grassland
# Levels: Desert < Grassland < Forest


#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$value, ktdat$group_label, method = "bh")
pt

pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )

sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a         a 
# 3 Grassland      b          b
# 2    Forest      b          b
sig <- sig[ order(sig$Group), ]

levels(sig$Group) # "Desert"    "Grassland" "Forest"   





####
#### Assess CPP in selected compounds
####

dat <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Liu_2020_mice_soil_PRJNA542998-Hostremoval.rds")

dat <- filter(dat, group_label %in% c("Desert","Grassland","Forest"))


## c) Glucose
sel.cpd <- which(df.comp$name == "D-Glucose")
this_var <- "Glucose"

df.comp[sel.cpd, ]
#             id    abbrev      name    form 
# 27    cpd00027     glc-D D-Glucose C6H12O6 
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6   

sel <- which(dat$cpd_id == "cpd00027")

head(dat[sel, ])

temp_dat <- dat[sel, ]

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.3","~x ","~10^{-13}") #

p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(0.8, 1.8, 2.8), y = c(0.85, 0.6, 0.75), label = c("a","c","b"), size = 3.75)+
  
  labs(x = NULL, y = "Glucose CPP(%)") +
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-c-Glucose-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")

ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 59.38, df = 2, p-value = 1.276e-13
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.3","~x ","~10^{-13}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a        a  
# 3 Grassland      c          c
# 2    Forest      b         b 
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  






## d) Cellulose
sel.cpd <- which(df.comp$name == "Cellulose")
this_var <- "Cellulose"

df.comp[sel.cpd, ]
#             id    abbrev      name      form  OC_ratio HC_ratio NC_ratio 
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2 0.8333333 1.666667        0

sel <- which(dat$cpd_id == "cpd11746")

head(dat[sel, ])
temp_dat <- dat[sel, ]

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.2","~x ","~10^{-9}") #

p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(1.2, 2.2, 3.2), y = c(0.0026, 0.002, 0.0026), label = c("a","b","a"), size = 3.75)+
  
  labs(x = NULL, y = "Cellulose CPP(%)") +
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 12.5,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-d-Cellulose-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")

ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 41.102, df = 2, p-value = 1.188e-09
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.2","~x ","~10^{-9}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
#      Group Letter MonoLetter
# 1    Desert      a         a 
# 3 Grassland      b          b
# 2    Forest      a         a 
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  







## e) CO2 - "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2

sel <- which(dat$cpd_id == "cpd00011")


head(dat[sel, ])
temp_dat <- dat[sel, ]

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 3.8","~x ","~10^{-14}") #

p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  labs(x = NULL, y = "CO2 CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "bottom", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(1.2, 2.2, 3.2), y = c(0.59, 0.73, 0.73), label = c("a","b","b"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(e)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-e-CO2-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 61.8, df = 2, p-value = 3.805e-14
ktresult <- paste0("~Kruskal-Wallis:","~P == 3.8","~x ","~10^{-14}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a         a 
# 3 Grassland      b          b
# 2    Forest      b          b
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  




## f) will be heatmap



# ## g) O2 - "Oxygen"
# sel.cpd <- which(df.comp$name == "O2")
# this_var <- "O2"
# 
# df.comp[sel.cpd, ]
# #         id abbrev name form
# # 7 cpd00007     o2   O2   O2   
# 
# sel <- which(dat$cpd_id == "cpd00007")
# 
# head(dat[sel, ])
# temp_dat <- dat[sel, ]
# 
# 
# p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
#   #theme_bw()+
#   theme_classic()+
#   geom_boxplot(outlier.shape = NA)+
#   #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
#   geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
#   scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
#   
#   labs(x = NULL, y = "O2 CPP(%)") +
#   
#   theme(
#     #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
#     #legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill = "transparent"),
#     legend.key.size = unit(0.9,"line"),
#     legend.title = element_text(size = rel(0.9)),
#     legend.text = element_text(size = rel(0.8)) ,
#     axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
#   )
# 
# p
# 
# 
# grid.text(label = "(g)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
# dev.print(tiff, file = paste0(workdir,"/plots/","CPP-g-O2-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")
# 

## instead compare H2O ??

## g) H2O - "Water"
sel.cpd <- which(df.comp$name == "H2O")
this_var <- "H2O"

df.comp[sel.cpd, ]
#         id abbrev name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# 1 cpd00001    h2o  H2O  H2O 

sel <- which(dat$cpd_id == "cpd00001")

head(dat[sel, ])
temp_dat <- dat[sel, ]

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.6","~x ","~10^{-15}") #


p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  labs(x = NULL, y = "H2O CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(0.8, 1.8, 2.8), y = c(5.6, 4.7, 5.2), label = c("a","c","b"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p


grid.text(label = "(g)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-g-H2O-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")

ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 68.167, df = 2, p-value = 1.576e-15
ktresult <- paste0("~Kruskal-Wallis:","~P == 1.6","~x ","~10^{-15}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a        a  
# 3 Grassland      c          c
# 2    Forest      b         b 
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  




## h) 
# adenylate energy charge (AEC) indicates the energetic status of soil microorganisms
# the energy status of soilmicroorganisms was evaluated by determining AEC defined as: 
# AEC = (ATP + 0.5 × ADP) / (ATP + ADP + AMP)

sel.cpd <- which(df.comp$name == "ATP")
df.comp[sel.cpd, ]
#.        id abbrev name          form OC_ratio HC_ratio NC_ratio 
# 2 cpd00002    atp  ATP C10H13N5O13P3      1.3      1.3      0.5  
sel <- which(dat$cpd_id == "cpd00002")
head(dat[sel, ])
vals <- list()
vals[["ATP"]] <- dat[sel, ]
vals[["ATP"]]$sample

sel.cpd <- which(df.comp$name == "ADP")
df.comp[sel.cpd, ]
#         id abbrev name          form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio
# 8 cpd00008    adp  ADP C10H13N5O10P2        1      1.3      0.5 
sel <- which(dat$cpd_id == "cpd00008")
head(dat[sel, ])
vals[["ADP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["ADP"]]$sample ) # TRUE

sel.cpd <- which(df.comp$name == "AMP")
df.comp[sel.cpd, ]
#          id abbrev name        form OC_ratio HC_ratio NC_ratio
# 18 cpd00018    amp  AMP C10H12N5O7P      0.7      1.2      0.5 
sel <- which(dat$cpd_id == "cpd00018")
head(dat[sel, ])
vals[["AMP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["AMP"]]$sample ) # TRUE

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
AMP <- vals[["AMP"]]$cpd_rel_abun

AEC <- (ATP + 0.5*ADP) / (ATP + ADP + AMP)

temp <- cbind(dat[sel, ],data.frame(AEC=AEC))
head(temp)

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 4.7","~x ","~10^{-13}") #

p <- ggplot(data=temp, aes(x=group_label, y=AEC)) +
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  labs(x = NULL, y = "AEC, from CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(0.8, 1.8, 2.8), y = c(0.695, 0.675, 0.675), label = c("a","b","b"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(h)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-h-AEC-ratio-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp
# Kruskal-Wallis test
kt <- kruskal.test( AEC ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  AEC by group_label
# Kruskal-Wallis chi-squared = 56.789, df = 2, p-value = 4.66e-13
ktresult <- paste0("~Kruskal-Wallis:","~P == 4.7","~x ","~10^{-13}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$AEC, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a         a 
# 3 Grassland      b          b
# 2    Forest      b          b
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  



## i) ATP / ADP

# use vectors for ATP and ADP from above

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
#AMP <- vals[["AMP"]]$cpd_rel_abun

ATP_ADP_ratio <- ATP/ADP

temp <- cbind(dat[sel, ],data.frame(ATP_ADP_ratio=ATP_ADP_ratio))

# Include Kruskal-Wallis & post-hoc Dunn test result, per below
ktresult <- paste0("~Kruskal-Wallis:","~P == 6.1","~x ","~10^{-13}") #

p <- ggplot(data=temp, aes(x=group_label, y=ATP_ADP_ratio)) +
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
  
  labs(x = NULL, y = "ATP/ADP, from CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "bottom", label = ktresult, parse = TRUE, size = 3.25 , lineheight = 0.85 )+
  annotate("text", x = c(0.8, 1.8, 2.8), y = c(1.55, 1.71, 1.71), label = c("a","b","b"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(i)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-i-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp
# Kruskal-Wallis test
kt <- kruskal.test( ATP_ADP_ratio ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  ATP_ADP_ratio by group_label
# Kruskal-Wallis chi-squared = 56.236, df = 2, p-value = 6.143e-13
ktresult <- paste0("~Kruskal-Wallis:","~P == 6.1","~x ","~10^{-13}") #

## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$ATP_ADP_ratio, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "Desert"    "Forest"    "Grassland"
sig$Group <- factor( sig$Group, levels=c("Desert", "Grassland", "Forest"),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1    Desert      a         a 
# 3 Grassland      b          b
# 2    Forest      b          b
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "Desert"    "Grassland" "Forest"  



## f) HEATMAP

dim(dat) # 709063      6
head(dat)
# cpd_id     sample cpd_rel_abun log10_abun            group group_label
# 1 cpd25681 SRR9276174 2.633358e-05 -4.5794901 F_first_sampling      Forest
# 2 cpd02597 SRR9276174 2.214690e-04 -3.6546871 F_first_sampling      Forest
# 3 cpd24620 SRR9276174 2.633358e-05 -4.5794901 F_first_sampling      Forest
# 4 cpd00001 SRR9276174 4.882249e+00  0.6886199 F_first_sampling      Forest
# 5 cpd01501 SRR9276174 2.074767e-04 -3.6830307 F_first_sampling      Forest
# 6 cpd00851 SRR9276174 1.344470e-03 -2.8714489 F_first_sampling      Forest

# p<- ggplot(dat, aes(x = sample, y = cpd_id, fill = log10_abun)) + # ggplot(long_df, aes(x = column_name, y = row_id, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_minimal()
# p


#library(pheatmap)

dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ sample', value.var = "log10_abun")
#dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ abbrev', value.var = "cpd_rel_abun")
row.names(dat.wide) <- dat.wide[ ,1]
dat.wide <- dat.wide[ ,-1]

# delete compounds with minimal variation
#sel <- which(dat.wide == 0) # 

# calculate row standard deviation?
row_sd <- apply(dat.wide, 1, sd, na.rm = TRUE)
hist(row_sd); summary(row_sd)
# cpd_rel_abun
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1358  0.3501  0.6213  1.0721  2.9434 


dim(dat.wide) # 7967   89
length(row_sd) # 7967
#low_sd_rows <- df[row_sd < 0.5, ]

quantile(row_sd, probs = 0.5)
# cpd_rel_abun
# 50% 
# 0.3501007 
sel <- which(row_sd > quantile(row_sd, probs = 0.5) ) # n = 3982 compounds

dat.wide.select <- dat.wide[sel, ]

colnames(dat.wide.select)
identical( colnames(dat.wide.select) , rownames(phy_in@sam_data) ) # TRUE

# col.df <- data.frame(sample=colnames(dat.wide.select), group_label=phy_in@sam_data$group_label)
# row.names(col.df) <- col.df$sample

col.df <- data.frame(Source=phy_in@sam_data$group_label)
row.names(col.df) <- colnames(dat.wide.select)

cols.ecosystem_type
# Desert      Grassland         Forest 
# "gold" "springgreen4"      "purple4" 

ann_colors = list(
  Source = c(Desert = "gold", Grassland = "springgreen4", Forest = "purple4" )
)

#pheatmap(as.matrix(dat.wide))
#pheatmap(as.matrix(dat.wide), scale = "row", show_rownames = FALSE)

#pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 6)
pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 10)

pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 10, show_colnames = FALSE, annotation_col = col.df , annotation_colors = ann_colors )


grid.text(label = "(f)", x = unit(0.033, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=18, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-f-HEATMAP-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 18, height = 22, units = "cm", res = 500, compression = "lzw",type="cairo")


#-------------------------



##########################
##########################
##########################
##########################

## Supplementary validation dataset

#### Chuckran soil with glucose addition
#    prepare for data download from SRA
#-------------------------


sradat <- read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/metadata_table1_copy_chuckran_soil_glucose.xlsx", sheet = 1, range = "A1:P30")
sradat <- as.data.frame(sradat)
str(sradat)
# 'data.frame':	29 obs. of  16 variables:
#  $ IMG_id                   : num  3.3e+09 3.3e+09 3.3e+09 3.3e+09 3.3e+09 ...
# $ Sample_name              : chr  "C0D1" "C0D2" "C0D3" "C0D4" ...
# $ Time_h                   : num  0 0 0 0 0 0 0 0 24 24 ...
# $ JGI_analysis_project_type: chr  "Metagenome" "Metagenome" "Metagenome" "Metagenome" ...
# $ NCBI_BioProject_accession: chr  "PRJNA539715" "PRJNA539712" "PRJNA539720" "PRJNA539713" ...
# $ NCBI_BioSample_accession : chr  "SAMN11528526" "SAMN11533409" "SAMN11533240" "SAMN11532952" ...
# $ SRA_accession            : chr  "SRR9032300" "SRR9032202" "SRR9032615" "SRR9032258" ...
# $ No_reads                 : num  1.23e+08 1.11e+08 1.26e+08 1.29e+08 9.90e+07 ...
# $ Assembled_genome_size_bp : num  1.36e+09 1.15e+09 1.65e+09 1.47e+09 3.03e+08 ...
# $ No_of_genes              : num  3466757 2964388 4265445 3759875 718971 ...
# $ No_of_scaffolds          : num  2855912 2455617 3536653 3098631 631914 ...
# $ N50_bp                   : num  851672 739392 1067790 927590 217331 ...
# $ GC_content               : num  63.9 64.1 64.4 64.1 60.4 57.4 58.7 59.2 64 64.3 ...
# $ COG_database             : num  53.7 53.9 53.9 53.8 45.8 37.6 41 46.7 52.6 54.4 ...
# $ Pfam_database            : num  51.5 51.6 51.7 51.6 43.2 35.3 38.6 43.6 51.5 51.8 ...
# $ KEGG_database            : num  24.9 25 25.3 25 17.8 14.8 15.9 18 24.4 25.6 ...

# only keep metagenomes

sel <- which(sradat$JGI_analysis_project_type == "Metagenome") # 13

table( sradat$Time_h[sel] )
# 0  8 24 48 
# 4  3  3  3 



sradat.select <- sradat[sel, ]

dim(sradat.select) # 13 16

sradat.select$Run <- sradat.select$SRA_accession


temp <- sradat.select

saveRDS(object = sradat.select, file = "sradat.select.chuckran-soil-glucose.rds")


runlist <- read.csv( "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/chuckran-soil-glucose-runlist.txt", header = FALSE)
runlist <- runlist[,1]


identical(runlist, sradat.select$Run) # TRUE


#-------------------------

#### Chuckran soil with glucose addition - read in superfocus - fxn potential outputs
#-------------------------

sradat.select <- readRDS("sradat.select.chuckran-soil-glucose.rds")

sampid <- sradat.select$Run # 13

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032199"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032202"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032258"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032259"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032267"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032300"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR9032199 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032199" 
# SRR9032202 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032202" 
# SRR9032258 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032258" 
# SRR9032259 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032259" 
# SRR9032267 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032267" 
# SRR9032300 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/3_fxn_superfocus_copy/superfocus_out_SRR9032300" 

length(results_dirs) # 13

sel <- which(names(results_dirs) %in% sampid) # qty 13
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 13

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 13

# reset sampid to remove missing sample
sampid <- names(results_dirs)
identical(sampid, names(results_dirs)) # TRUE


# In this data one Run corresponds to a single Sample_ID !!!

# collate results into a long-format table

sfx.long <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,percent_abun=NA)

for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  # [1] "Subsystem.Level.1"                                                                         
  # [2] "Subsystem.Level.2"                                                                         
  # [3] "Subsystem.Level.3"                                                                         
  # [4] "Function"                                                                                  
  # [5] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq"  
  # [6] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq.."
  
  
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [7] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [8] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [9] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [10] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  # [11] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [12] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.fastq..$", x = names(tab)) # for single (unpaired) reads
  if (length(sel.col.percent)>1) {
    tab$percent_abun <- apply(X = tab[ ,sel.col.percent], MARGIN = 1, FUN = mean )
  } else {
    tab$percent_abun <- tab[ ,sel.col.percent]
  }
  
  # sum(tab$percent_abun) # 100
  # mean(tab$percent_abun) # 0.004338583
  
  names(sfx.long) # "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"    "percent_abun"
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # ...
  # [13] "sampid"
  # [14] "percent_abun"
  
  tab <- tab[ ,c("sampid","Subsystem.Level.1","Subsystem.Level.2","Subsystem.Level.3","Function","percent_abun")]
  names(tab) <- names(sfx.long)
  
  sfx.long <- rbind(sfx.long,tab)
  
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}


head(sfx.long)
# remove empty 1st row
sfx.long <- sfx.long[-1, ]
dim(sfx.long) # 391483      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2           subsys_L3
# 2 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# 3 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# 4 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# 5 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# 6 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# 7 SRR9032199 Amino Acids and Derivatives         - Amino acid racemase
# fxn percent_abun
# 2                                         2-methylaconitate_cis-trans_isomerase 5.602917e-04
# 3                                                   2-methylaconitate_isomerase 7.204053e-06
# 4 2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117) 1.248349e-04
# 5                       2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79) 1.130401e-04
# 6                                       4-hydroxyproline_epimerase_(EC_5.1.1.8) 2.259036e-03
# 7                                    4-oxalomesaconate_tautomerase_(EC_5.3.2.8) 5.226470e-05

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 36025    14

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 36025
length(unique(full_fxn_names)) # 36025

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_isomerase" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-oxalomesaconate_tautomerase_(EC_5.3.2.8)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR9032199"   "SRR9032202"   "SRR9032258"   "SRR9032259"   "SRR9032267"   "SRR9032300"  
# [8] "SRR9032509"   "SRR9032510"   "SRR9032615"   "SRR9032617"   "SRR9032694"   "SRR9032715"   "SRR9032716" 

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" "SRR9032267" "SRR9032300"

length(sampid) # 13

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # TRUE
#identical(sort(names(sfx.wide)), sort(as.character(sampid))) #

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 36025     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 194
length(unique(tax.fxn$subsys_L3)) # 1270
length(unique(tax.fxn$fxn)) # 17012


# # # #

## gather Function count data??
sfx.long.count <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,count_abun=NA)
length(sampid) # 13
for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  tab$sampid <- this_samp
  names(tab)
  tab <- tab[,c(7,1,2,3,4,5)] # this time capture 'count' data
  names(tab) <- names(sfx.long.count)
  sfx.long.count <- rbind(sfx.long.count,tab)
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}
head(sfx.long.count)
# remove empty 1st row
sfx.long.count <- sfx.long.count[-1, ]
sum(sfx.long.count$count_abun) # 258744401 = 258,744,401
sfx.long.count$full_fxn_tax <- paste0(sfx.long.count$subsys_L1,"___", sfx.long.count$subsys_L2,"___", sfx.long.count$subsys_L3,"___", sfx.long.count$fxn)
head(sfx.long.count)
sfx.wide.count <- dcast(sfx.long.count, formula = full_fxn_tax ~ sampleID, value.var = "count_abun")
dim(sfx.wide.count) # 36025    14
sel.na <- which(is.na(sfx.wide.count),arr.ind = TRUE)
sfx.wide.count[sel.na] <- 0
sum(colSums(sfx.wide.count[,-1])) # 258744401
hist(colSums(sfx.wide.count[,-1]))
mean(colSums(sfx.wide.count[,-1])) # 19903415
sd(colSums(sfx.wide.count[,-1])) # 3851888

summary(colSums(sfx.wide.count[,-1]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 15255213 17182086 18674078 19903415 22024294 28954895 
length(unique(sfx.long.count$subsys_L1)) # 35

fxn_sum_counts <- colSums(sfx.wide.count[,-1])

# # # #


#-------------------------

#### Chuckran soil with glucose addition - functions - get into Phyloseq object
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 36025     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 36025    13

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 36025 taxa and 13 samples ]
# tax_table()   Taxonomy Table:    [ 36025 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" "SRR9032267" "SRR9032300" "SRR9032509" "SRR9032510"
# [9] "SRR9032615" "SRR9032617" "SRR9032694" "SRR9032715" "SRR9032716"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

head(row.names(sradat.select))

samp <- sradat.select

dim(samp) # 13 17

head(row.names(samp)) # 

row.names(samp) <- samp$Run

identical(row.names(samp), sample_names(phy)) # FALSE
length(row.names(samp)) # 13
length(sample_names(phy)) # 13
sel <- which(row.names(samp) %in% sample_names(phy)) # 13

samp2 <- samp[sample_names(phy), ]

identical(row.names(samp2), names(fxn_sum_counts)) # TRUE

samp2$fxn_sum_counts <- fxn_sum_counts


SAMP <- sample_data(samp2)



### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 36025 taxa and 13 samples ]
# sample_data() Sample Data:       [ 13 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 36025 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase"
# fxn                                                                            
# fxn_1 "2-methylaconitate_cis-trans_isomerase"                                        
# fxn_2 "2-methylaconitate_isomerase"                                                  
# fxn_3 "2-methylcitrate_dehydratase_(2-methyl-trans-aconitate_forming)_(EC_4.2.1.117)"
# fxn_4 "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                      
# fxn_5 "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                      
# fxn_6 "4-oxalomesaconate_tautomerase_(EC_5.3.2.8)"   


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = phy, file = "phy-phyloseq-fxn-chuckran-soil-glucose.RDS")


head(phy@sam_data)

# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 13
NonZeroFxns

mean(NonZeroFxns) # 30114.08
sd(NonZeroFxns) # 722.9469


table(phy@sam_data$Time_h)
# 0  8 24 48 
# 4  3  3  3 

class(phy@sam_data$Time_h) # "numeric"

phy@sam_data$group_label <- factor(phy@sam_data$Time_h, levels = c(0,8,24,48), labels = c("0h","8h","24h","48h"),ordered=TRUE)


#-------------------------

#### Chuckran soil with glucose addition - COPY of R code to run CPP steps on HPC
#    1) build reaction search - get reactions and compounds
#    2) get cpd rel abun per sample
#    3) collate compounds for each sample
#-------------------------

# # # # # # # # # # # # #
# #
# # R script for cpp3d
# # - build reaction search in parallel - get_reactions & compounds
# # - get cpd rel abun per sample
# # - collate_compounds
# #
# # For study - Chuckran et al soil with glucose addition
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# message("\n# establish folders and input files")
# 
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis'")
# workdir <- "/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis/working"
# 
# message("\nthis_study <- '-chuckran-soil-glucose-pawsey'")
# this_study <- "-chuckran-soil-glucose-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-chuckran-soil-glucose.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-chuckran-soil-glucose.RDS")
# 
# 
# subsys.lut <- readRDS("subsys.lut.RDS")
# rxns.lut <- readRDS("rxns.lut.RDS")
# rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
# compounds.lut <- readRDS("compounds.lut.RDS")
# 
# 
# 
# message("\n### 1) build reaction search in parallel - get_reactions & compounds")
# message("\n# # # # # # # # # #")
# message("\ndf.tax <- as.data.frame(phy@tax_table)")
# df.tax <- as.data.frame(phy@tax_table)
# message("\nhead(row.names(df.tax))")
# head(row.names(df.tax))
# message("\ndim(df.tax)")
# dim(df.tax)
# 
# 
# get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
#   
#   rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#   rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#   
#   sub1 <- df.tax$subsys_L1[i]
#   sub2 <- df.tax$subsys_L2[i]
#   sub3 <- df.tax$subsys_L3[i]
#   
#   fxn.temp <- df.tax$fxn[i]
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   
#   # store results corresponding to each Superfocus row
#   fxn.list <- list()
#   fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
#   
#   # check for multiple functions/reactions?
#   flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
#   flag2 <- grepl(pattern = "_@_", x = fxn.temp)
#   if (!any(flag1,flag2)==TRUE) {
#     # no multiples
#     fxns <- fxn.temp
#   } else if (flag1==TRUE) {
#     fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
#   } else {
#     fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
#   }
#   # remove underscores
#   ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
#   
#   # process each fxn & store attributes
#   df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
#   
#   # Identify '/' separators with no '_'  ??
#   
#   for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
#     #f<-1
#     #f<-2
#     f.in <- fxns[f]
#     
#     # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
#     full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
#     full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
#     
#     ## data cleaning
#     
#     # trim off '_#' and '_##' tags
#     trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
#     if (!is.na(trim_nchar) & length(trim_nchar)==1) {
#       f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
#     }
#     
#     # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
#     f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
#     
#     #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#     #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#     
#     full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
#     full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
#     
#     sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
#     
#     ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
#     if (length(sel.rx)>=1) {
#       df.fxns$matching_method[f] <- "Exact hierachy match"
#       df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#       
#     } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
#       # search by EC id if present
#       
#       f.in <- fxns[f] # this goes back to string with brackets for EC
#       ## LOOK FOR MULTIPLE ECs ?
#       
#       how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
#       
#       ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
#       #class(ECs)
#       ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
#       ECs.collapse <- paste0(ECs, collapse = "|")
#       
#       sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
#       
#       if (length(how_many_ECs)==0 | length(ECs)==0) {
#         # there was a glitch, database typo, or some error in identifying the EC number
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       } else if (length(sel.rx)>=1) {
#         # combined EC hits identified
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
#         # treat EC hits individually
#         sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
#         # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
#         sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else {
#         # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       # END EC matching
#       
#       
#     } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
#       ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
#       # 1. 'reactions' table by name: rxns.lut$name
#       # 2. 'reactions' table by aliases: rxns.lut$aliases
#       # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
#       # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
#       
#       if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
#         # 1a - exact match - rxns.lut$name
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
#         #rxns.lut$name[sel.rx]
#         df.fxns$matching_method[f] <- "Matched Reactions name"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
#         # 2a - exact match - rxns.lut$aliases
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
#         #rxns.lut$aliases[sel.rx]
#         #rxns.lut$name[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Reactions aliases"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
#         # 3a - exact match - subsys.lut$Role
#         sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
#         #subsys.lut$Role[sel.rx]
#         #subsys.lut$Reaction[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Subsytem role"
#         df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
#         # 4a - exact match - rxn_pathways.lut$External_rxn_name
#         sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
#         
#         df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
#         df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
#         
#         
#       } else {
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       
#       ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
#       
#       
#     } # END function - reaction search
#     
#     #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
#     
#     #print(paste0("completed fxn ", f))
#     
#     
#     ## now investigate these reactions ...
#     # Reactions lookup table: 
#     # - "equation": Definition of reaction expressed using compound IDs and after protonation
#     # Compounds lookup table:
#     # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
#     #df.fxns
#     
#     
#     #if (df.fxns$matching_method == "No match found") {
#     if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
#       
#       df.Rxns <- NA
#       df.Compounds <- NA
#       
#     } else { # reaction(s) were identified
#       
#       # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
#       f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
#       
#       df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
#                             rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
#       
#       for (r in 1:dim(df.Rxns)[1]) {
#         #r<-1
#         #this_rxn <- "rxn00004"
#         this_rxn <- df.Rxns$rxn_id[r]
#         sel <- which(rxns.lut$id == this_rxn)
#         ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
#         ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
#         ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
#         
#         # extract compound info
#         
#         #df.Rxns$rxn_eqn[r]
#         #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
#         #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
#         
#         ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
#         # 5 23 43 61
#         # 6 25 46 65 83
#         
#         ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
#         # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
#         
#         if (length(compds)>=1) {
#           
#           df.Rxns$compds[r] <- paste0(compds, collapse = ";")
#           
#           ## get compound coefficients?
#           start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
#           end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
#           ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
#           
#           df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
#           
#           # get formulas of compounds
#           
#           formx <-filter(compounds.lut, id %in% compds )
#           row.names(formx) <- formx$id
#           ( formx.char <- formx[compds, ]$formula )
#           # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
#           # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
#           # "C7H7O7" "H2O"    "C7H5O6"
#           df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
#           
#           ( compd.names <- formx[compds, ]$name )
#           # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
#           
#           temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
#                                           cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
#                                           #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
#           )
#           
#         } else {
#           # No specified reaction equation or chemical formula info
#           df.Rxns$compds[r] <- NA
#           df.Rxns$compd_coef[r] <- NA
#           df.Rxns$chem_formx[r] <- NA
#           
#           temp.df.Compounds <- NA
#           
#         }
#         
#         if (r==1) { df.Compounds <- temp.df.Compounds }
#         
#         if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
#         
#         # clean up - if there are additional reactions?
#         temp.df.Compounds <- NA
#         
#       } # END loop for r - rxn_id's per f/f.in
#       
#     } # END else loop when reactions identified
#     
#     # store results corresponding to each sub-reaction of each Superfocus row
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
#     
#     
#   } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
#   
#   saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
#   
# } # END function to be run in parallel for each superfocus row
# 
# 
# # # # # # # # # # # # # # # # # # #
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
#   get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# message("\n(num_results_files <- dim(df.tax)[1])")
# (num_results_files <- dim(df.tax)[1])
# 
# # assemble all compound data outputs
# # start with blank row
# 
# df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )
# 
# for (i in 1:num_results_files) {
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
#   
#   f_no <- length( temp[[1]][["compounds"]] )
#   
#   for (f in 1:f_no) {
#     #f<-2
#     # only add non-NA results
#     if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
#       
#       df.temp <- temp[[1]][["compounds"]][[f]]
#       ok <- complete.cases(df.temp)
#       df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
#       df.out <- rbind(df.out,df.temp)
#     }
#   }
#   print(paste0("added df ",i," of ",num_results_files ))
#   
# }
# 
# 
# message("\nstr(df.out)")
# str(df.out)
# 
# 
# saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))
# 
# # remove NA first row
# message("\nhead(df.out)")
# head(df.out)
# 
# df.out <- df.out[-1, ]
# 
# message("\ndim(df.out)")
# dim(df.out)
# 
# 
# message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")
# 
# df.out$cpd_molar_prop_norm <- NA
# 
# message("\nlength(unique(df.out$superfocus_fxn))")
# length(unique(df.out$superfocus_fxn))
# 
# message("\nphy")
# phy
# 
# message("\n% of functions represented - with compound information")
# 100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))
# 
# 
# fxns_found <- unique(df.out$superfocus_fxn)
# 
# for (k in 1:length(fxns_found)) {
#   #k<-1
#   this_fxn <- fxns_found[k]
#   sel <- which(df.out$superfocus_fxn == this_fxn)
#   
#   sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
#   # calculate 
#   
#   df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
#   
#   print(paste0("completed ",k))
#   
# }
# 
# message("\nsum(df.out$cpd_molar_prop_norm)")
# sum(df.out$cpd_molar_prop_norm)
# 
# message("\nsample_sums(phy)")
# sample_sums(phy)
# 
# message("\ngetwd()")
# getwd()
# 
# saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))
# 
# 
# 
# message("\n### 2) get cpd rel abun per sample")
# message("\n# # # # # # # # # #")
# 
# 
# df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
# message("\ndim(df.OTU)")
# dim(df.OTU)
# 
# 
# get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
#   
#   this_samp <- sample_names(phy_in)[i]
#   df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
#   
#   dat.cpd$sample <- this_samp
#   
#   dat.cpd$cpd_rel_abun_norm <- NA
#   
#   fxns_all <- row.names(df.OTU)
#   
#   for (k in 1:length(fxns_all)) {
#     #k<-1
#     this_fxn <- fxns_all[k]
#     sel <- which(dat.cpd$superfocus_fxn == this_fxn)
#     
#     if (length(sel)>=1) {
#       dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
#       
#     }
#   } # END rel abun values for all relevant functions added
#   
#   saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   
# } # END
# 
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
# head(dat)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   dat <- rbind(dat, temp)
#   print(paste0("completed ",i))
# }
# 
# 
# saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )
# 
# rm(temp)
# 
# message("\nstr(dat)")
# str(dat)
# 
# message("\nsum(dat$cpd_rel_abun_norm)")
# sum(dat$cpd_rel_abun_norm)
# 
# message("\naverage functional relative abundance per sample")
# message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# 
# message("\nnames(dat)")
# names(dat)
# 
# message("\nlength(unique(dat$cpd_id))")
# length(unique(dat$cpd_id))
# 
# 
# 
# 
# message("\n### 3) collate_compounds within each sample")
# message("\n# # # # # # # # # #")
# 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
# } # END
# 
# 
# 
# no_forks <- 4
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# message("\nstr(dat.cpd.collate)")
# str(dat.cpd.collate)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)")
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))
# 
# # END


#-------------------------

#### Chuckran soil with glucose addition - COPY of OUTOUTS from R code after running CPP steps on HPC
#-------------------------

# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "4.1"
# 
# $year
# [1] "2024"
# 
# $month
# [1] "06"
# 
# $day
# [1] "14"
# 
# $`svn rev`
# [1] "86737"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.4.1 (2024-06-14)"
# 
# $nickname
# [1] "Race for Your Life"
# 
# [1] ‘4.4.1’
# Loading required package: foreach
# Loading required package: iterators
# [1] ‘1.0.17’
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# [1] ‘1.1.4’
# [1] ‘1.5.2’
# [1] ‘1.46.0’
# 
# # establish folders and input files
# 
# workdir <- '/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis/working'
# 
# this_study <- '-chuckran-soil-glucose-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-chuckran-soil-glucose.RDS')
# 
# ### 1) build reaction search in parallel - get_reactions & compounds
# 
# # # # # # # # # # #
# 
# df.tax <- as.data.frame(phy@tax_table)
# 
# head(row.names(df.tax))
# [1] "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"
# 
# dim(df.tax)
# [1] 36025     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# ...
# 
# 
# 
# [[36023]]
# NULL
# 
# [[36024]]
# NULL
# 
# [[36025]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 36025
# [1] "added df 1 of 36025"
# [1] "added df 2 of 36025"
# [1] "added df 3 of 36025"
# ...
# 
# 
# [1] "added df 36023 of 36025"
# [1] "added df 36024 of 36025"
# [1] "added df 36025 of 36025"
# 
# str(df.out)
# 'data.frame':	1295095 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id        : chr  NA "rxn25278" "rxn25278" "rxn25279" ...
# $ cpd_id        : chr  NA "cpd25681" "cpd02597" "cpd24620" ...
# $ cpd_name      : chr  NA "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form      : chr  NA "C7H5O6" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f
# 1           <NA> NA
# 2          fxn_2  1
# 3          fxn_2  1
# 4          fxn_3  1
# 5          fxn_3  1
# 6          fxn_3  1
# f__in
# 1                                                                          <NA>
#   2                                                   2-methylaconitate isomerase
# 3                                                   2-methylaconitate isomerase
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 6 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn25278 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 3 rxn25278 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 4 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 5 rxn25279 cpd00001                                          H2O      H2O
# 6 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 1295094       8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 19178
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 36025 taxa and 13 samples ]
# sample_data() Sample Data:       [ 13 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 36025 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 53.23525
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 19175"
# [1] "completed 19176"
# [1] "completed 19177"
# [1] "completed 19178"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 19178
# 
# sample_sums(phy)
# SRR9032199 SRR9032202 SRR9032258 SRR9032259 SRR9032267 SRR9032300 SRR9032509 
# 100        100        100        100        100        100        100 
# SRR9032510 SRR9032615 SRR9032617 SRR9032694 SRR9032715 SRR9032716 
# 100        100        100        100        100        100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/chuckran_soil_glucose/cpp_analysis"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 36025    13
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# ...
# 
# 
# 
# [[12]]
# NULL
# 
# [[13]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f
# 2          fxn_2 1
# 3          fxn_2 1
# 4          fxn_3 1
# 5          fxn_3 1
# 6          fxn_3 1
# 7          fxn_4 1
# f__in
# 2                                                   2-methylaconitate isomerase
# 3                                                   2-methylaconitate isomerase
# 4 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 5 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 6 2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)
# 7                       2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn25278 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 3 rxn25278 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 4 rxn25279 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 5 rxn25279 cpd00001                                          H2O      H2O
# 6 rxn25279 cpd25681                     2-methyl-trans-aconitate   C7H5O6
# 7 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# cpd_molar_prop cpd_molar_prop_norm     sample cpd_rel_abun_norm
# 2              1          0.50000000 SRR9032199      3.602027e-06
# 3              1          0.50000000 SRR9032199      3.602027e-06
# 4              1          0.33333333 SRR9032199      4.161165e-05
# 5              1          0.33333333 SRR9032199      4.161165e-05
# 6              1          0.33333333 SRR9032199      4.161165e-05
# 7              1          0.05555556 SRR9032199      6.280004e-06
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# [1] "completed 5"
# [1] "completed 6"
# [1] "completed 7"
# [1] "completed 8"
# [1] "completed 9"
# [1] "completed 10"
# [1] "completed 11"
# [1] "completed 12"
# [1] "completed 13"
# 
# str(dat)
# 'data.frame':	16836222 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylaconitate isomerase" "2-methylaconitate isomerase" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" "2-methylcitrate dehydratase (2-methyl-trans-aconitate forming) (EC 4.2.1.117)" ...
# $ rxn_id             : chr  "rxn25278" "rxn25278" "rxn25279" "rxn25279" ...
# $ cpd_id             : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ cpd_name           : chr  "2-methyl-trans-aconitate" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" "H2O" ...
# $ cpd_form           : chr  "C7H5O6" "C7H5O6" "C7H7O7" "H2O" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.333 0.333 0.333 ...
# $ sample             : chr  "SRR9032199" "SRR9032199" "SRR9032199" "SRR9032199" ...
# $ cpd_rel_abun_norm  : num  3.60e-06 3.60e-06 4.16e-05 4.16e-05 4.16e-05 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 867.1734
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 66.70565
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 8486
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# [[4]]
# NULL
# 
# [[5]]
# NULL
# 
# [[6]]
# NULL
# 
# [[7]]
# NULL
# 
# [[8]]
# NULL
# 
# [[9]]
# NULL
# 
# [[10]]
# NULL
# 
# [[11]]
# NULL
# 
# [[12]]
# NULL
# 
# [[13]]
# NULL
# 
# 
# ## assemble results
# cpd_id     sample cpd_rel_abun
# 1 cpd25681 SRR9032199 0.0004249709
# 2 cpd02597 SRR9032199 0.0191206785
# 3 cpd24620 SRR9032199 0.0005762276
# 4 cpd00001 SRR9032199 5.0011883693
# 5 cpd01501 SRR9032199 0.0131099817
# 6 cpd00851 SRR9032199 0.0108064371
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# [1] "completed 5"
# [1] "completed 6"
# [1] "completed 7"
# [1] "completed 8"
# [1] "completed 9"
# [1] "completed 10"
# [1] "completed 11"
# [1] "completed 12"
# [1] "completed 13"
# 
# str(dat.cpd.collate)
# 'data.frame':	110318 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9032199" "SRR9032199" "SRR9032199" "SRR9032199" ...
# $ cpd_rel_abun: num  0.000425 0.019121 0.000576 5.001188 0.01311 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 867.1734
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 66.70565
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries



#-------------------------

#### Chuckran soil with glucose addition - continue CPP analysis
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-chuckran-soil-glucose.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/chuckran-glucose-soil/cpp_analysis/dat.cpd.collate-all-samps-cpp3d--chuckran-soil-glucose-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	110318 obs. of  3 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9032199" "SRR9032199" "SRR9032199" "SRR9032199" ...
# $ cpd_rel_abun: num  0.000425 0.019121 0.000576 5.001188 0.01311 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000013 0.000167 0.007861 0.001501 8.012038 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -4.8752 -3.7767    -Inf -2.8236  0.9037 


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.8410 -4.8752 -3.7767 -4.0454 -2.8236  0.9037 


# make group variable from sample name

dat.cpd.collate$group <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 36025 taxa and 13 samples ]
# sample_data() Sample Data:       [ 13 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 36025 taxa by 4 taxonomic ranks ]

samp <- as(phy@sam_data, "data.frame")

unique(samp$Time_h)
#   8  0 24 48


for (i in 1:length(sample_names(phy))) {
  #for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  #this_samp <- samp$Run[i]
  sel.phy <- which(phy@sam_data$Run == this_samp)
  sel.dat <- which(dat.cpd.collate$sample == this_samp)
  
  dat.cpd.collate$group[sel.dat] <- phy@sam_data$Time_h[sel.phy]
  
  #dat.cpd.collate$group[sel.dat] <- as.character( samp$group_new[i] )
  
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group) #
#  8  0 24 48
class(dat.cpd.collate$group) # "numeric"

dat.cpd.collate$group_label <- factor(dat.cpd.collate$group, 
                                      levels = c( 0 , 8 , 24 , 48 ),
                                      labels = c("0h" , "8h" , "24h" , "48h" ),
                                      ordered = TRUE)

head(dat.cpd.collate)
# cpd_id     sample cpd_rel_abun log10_abun group group_label
# 1 cpd25681 SRR9032199 0.0004249709 -3.3716409     8          8h
# 2 cpd02597 SRR9032199 0.0191206785 -1.7184967     8          8h
# 3 cpd24620 SRR9032199 0.0005762276 -3.2394059     8          8h
# 4 cpd00001 SRR9032199 5.0011883693  0.6990732     8          8h
# 5 cpd01501 SRR9032199 0.0131099817 -1.8823979     8          8h
# 6 cpd00851 SRR9032199 0.0108064371 -1.9663175     8          8h

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-chuckran-soil-glucose.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-chuckran-soil-glucose.rds")


str(dat.cpd.collate)
# 'data.frame':	110318 obs. of  6 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9032199" "SRR9032199" "SRR9032199" "SRR9032199" ...
# $ cpd_rel_abun: num  0.000425 0.019121 0.000576 5.001188 0.01311 ...
# $ log10_abun  : num  -3.372 -1.718 -3.239 0.699 -1.882 ...
# $ group       : num  8 8 8 8 8 8 8 8 8 8 ...
# $ group_label : Ord.factor w/ 4 levels "0h"<"8h"<"24h"<..: 2 2 2 2 2 2 2 2 2 2 ...


length( unique(dat.cpd.collate$cpd_id) ) # 8486
8486*13 # 110318


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
#     cpd_id     sample cpd_rel_abun log10_abun group group_label
# 1 cpd25681 SRR9032199 0.0004249709 -3.3716409     8          8h
# 2 cpd02597 SRR9032199 0.0191206785 -1.7184967     8          8h
# 3 cpd24620 SRR9032199 0.0005762276 -3.2394059     8          8h
# 4 cpd00001 SRR9032199 5.0011883693  0.6990732     8          8h
# 5 cpd01501 SRR9032199 0.0131099817 -1.8823979     8          8h
# 6 cpd00851 SRR9032199 0.0108064371 -1.9663175     8          8h

dim(data_in) # 110318      6

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 8304.308
sd(no_compounds) #  23.62338

mean(sample_sum_relabun) # 66.70565
sd(sample_sum_relabun) # 0.1088342

#length(unique(data_in$cpd_id)) # 
length(unique(data_in$cpd_id[ which(data_in$cpd_rel_abun > 0) ])) # 8486


#-------------------------

#### Chuckran soil with glucose addition - continue CPP analysis
#    CPP - get into phyloseq object
#    beta diversity
#    alpha diversity
#    response in selected compounds: glucose, cellulose, CO2, O2, AEC, ATP/ADP
#    heatmap of scaled CPP
#-------------------------

phy.fxn <- readRDS("phy-phyloseq-fxn-chuckran-soil-glucose.RDS")
phy.fxn
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 36025 taxa and 13 samples ]
# sample_data() Sample Data:       [ 13 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 36025 taxa by 4 taxonomic ranks ]

table(phy.fxn@sam_data$Time_h)
# 0  8 24 48 
# 4  3  3  3

str(phy.fxn@sam_data)
# 'data.frame':	13 obs. of  18 variables:
#   Formal class 'sample_data' [package "phyloseq"] with 4 slots
# ..@ .Data    :List of 18
# .. ..$ : num  3.3e+09 3.3e+09 3.3e+09 3.3e+09 3.3e+09 ...
# .. ..$ : chr  "C8D2" "C0D2" "C0D4" "C8D1" ...
# .. ..$ : num  8 0 0 8 8 0 24 24 0 24 ...
# .. ..$ : chr  "Metagenome" "Metagenome" "Metagenome" "Metagenome" ...
# .. ..$ : chr  "PRJNA539711" "PRJNA539712" "PRJNA539713" "PRJNA539714" ...
# .. ..$ : chr  "SAMN11533337" "SAMN11533409" "SAMN11532952" "SAMN11532583" ...
# .. ..$ : chr  "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" ...
# .. ..$ : num  1.72e+08 1.11e+08 1.29e+08 1.32e+08 1.25e+08 ...
# .. ..$ : num  2.50e+09 1.15e+09 1.47e+09 1.84e+09 1.41e+09 ...
# .. ..$ : num  6199898 2964388 3759875 4756397 3629789 ...
# .. ..$ : num  5009710 2455617 3098631 3945810 2995887 ...
# .. ..$ : num  1429077 739392 927590 1187662 896872 ...
# .. ..$ : num  64.6 64.1 64.1 64.1 63.7 63.9 64 64.3 64.4 63.9 ...
# .. ..$ : num  53.4 53.9 53.8 52.6 53.2 53.7 52.6 54.4 53.9 52.9 ...
# .. ..$ : num  51.9 51.6 51.6 50.8 51.2 51.5 51.5 51.8 51.7 51.4 ...
# .. ..$ : num  24.7 25 25 24.4 24.7 24.9 24.4 25.6 25.3 24.7 ...
# .. ..$ : chr  "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" ...
# .. ..$ : num  23597826 15255213 17657472 18135466 17182086 ...
# ..@ names    : chr  "IMG_id" "Sample_name" "Time_h" "JGI_analysis_project_type" ...
# ..@ row.names: chr  "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" ...
# ..@ .S3Class : chr "data.frame"


samp <- as( phy.fxn@sam_data, "data.frame")
samp$group_label <- factor( samp$Time_h , levels = c(0,8,24,48), labels = c("0h","8h","24h","48h"), ordered = TRUE )
temp <- samp


dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-chuckran-soil-glucose.rds")

str(dat.cpd.collate)
# 'data.frame':	110318 obs. of  6 variables:
# $ cpd_id      : chr  "cpd25681" "cpd02597" "cpd24620" "cpd00001" ...
# $ sample      : chr  "SRR9032199" "SRR9032199" "SRR9032199" "SRR9032199" ...
# $ cpd_rel_abun: num  0.000425 0.019121 0.000576 5.001188 0.01311 ...
# $ log10_abun  : num  -3.372 -1.718 -3.239 0.699 -1.882 ...
# $ group       : num  8 8 8 8 8 8 8 8 8 8 ...
# $ group_label : Ord.factor w/ 4 levels "0h"<"8h"<"24h"<..: 2 2 2 2 2 2 2 2 2 2 ...

data_in <- dat.cpd.collate


str(data_in)
# 'data.frame':	110318 obs. of  6 variables:

length( unique(data_in$cpd_id) ) #  8486
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 8486
length( unique(data_in$sample) ) # 13



### get data into phyloseq object ...

head(data_in)
#     cpd_id     sample cpd_rel_abun log10_abun group group_label
# 1 cpd25681 SRR9032199 0.0004249709 -3.3716409     8          8h
# 2 cpd02597 SRR9032199 0.0191206785 -1.7184967     8          8h
# 3 cpd24620 SRR9032199 0.0005762276 -3.2394059     8          8h
# 4 cpd00001 SRR9032199 5.0011883693  0.6990732     8          8h
# 5 cpd01501 SRR9032199 0.0131099817 -1.8823979     8          8h
# 6 cpd00851 SRR9032199 0.0108064371 -1.9663175     8          8h


df.wide <- dcast(data_in, formula = sample + group_label ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#       sample group_label cpd00001 cpd00002 cpd00003 cpd00004  cpd00005  cpd00006  cpd00007 cpd00008
# 1 SRR9032199          8h 5.001188 2.329666 1.076612 1.047560 0.7396011 0.7422493 0.3227162 1.405330
# 2 SRR9032202          0h 4.978722 2.320058 1.065993 1.036932 0.7333483 0.7360414 0.3152727 1.395337
# 3 SRR9032258          0h 4.983385 2.324792 1.065785 1.036635 0.7346541 0.7372471 0.3173322 1.400866
# 4 SRR9032259          8h 4.993149 2.321264 1.070906 1.041768 0.7336205 0.7362626 0.3198508 1.396297
# 5 SRR9032267          8h 4.969843 2.315483 1.064068 1.034987 0.7308258 0.7334766 0.3165358 1.394570

# save group variable
samp <- df.wide[ ,1:2]
#samp <- df.wide[ ,1:3]
row.names(samp) <- samp$sample

# transpose
df.wide <- t(df.wide[ ,-2]) # minus 'group_label' column
#df.wide <- t(df.wide[ ,-c(2,3)]) #

head(df.wide)

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array"
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  8486   13
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix.
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)


# # convert Taxonomy table to matrix

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8486 taxa and 13 samples ]
# tax_table()   Taxonomy Table:    [ 8486 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# [1] "SRR9032199" "SRR9032202" "SRR9032258" "SRR9032259" "SRR9032267" "SRR9032300" "SRR9032509" "SRR9032510" "SRR9032615" "SRR9032617" "SRR9032694" "SRR9032715" "SRR9032716"

#identical(sample_names(phy.cpp), samp$sample) # TRUE
#identical(sample_names(phy.cpp), sradat.select2$Run) # TRUE
identical(sample_names(phy.cpp), temp$Run) # TRUE

#row.names(sradat.select2) <- sradat.select2$Run
identical( row.names(temp), sample_names(phy.cpp) ) # TRUE

# use earlier sample metadata from phy.fxn data object
#samp <- sradat.select2
samp <- temp

# row.names need to match sample_names() from phyloseq object
#row.names(samp) <- samp$sample
#identical(row.names(samp), samp$sample) # TRUE



### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8486 taxa and 13 samples ]
# sample_data() Sample Data:       [ 13 samples by 19 sample variables ]
# tax_table()   Taxonomy Table:    [ 8486 taxa by 1 taxonomic ranks ]

# check for 'taxa' (compounds) with zero data - because soil samples were excluded
min(taxa_sums(phy.cpp)) # 1.17471e-08
# # prune taxa that have zero sequence reads
# phy.cpp <- prune_taxa(taxa = taxa_sums(phy.cpp) > 0, x = phy.cpp)


phy.cpp@sam_data


saveRDS(object = phy.cpp, file = "phy.cpp-cleaned-chuckran-soil-glucose-v8b.RDS")

## remove sample "C0D3" due to outlying beta diversity ordination position

phy_in <- phy.cpp

sel <- which(phy_in@sam_data$Sample_name == "C0D3")

phy_in <- prune_samples(samples = sample_names(phy_in)[-sel], x = phy_in )
min(taxa_sums(phy_in)) # 0
# prune taxa that have zero sequence reads
phy_in <- prune_taxa(taxa = taxa_sums(phy_in) > 0, x = phy_in)

saveRDS(object = phy_in, file = "phy.cpp-cleaned-excludeC0D3-chuckran-soil-glucose-v8b.RDS")
phy.cpp <- readRDS("phy.cpp-cleaned-excludeC0D3-chuckran-soil-glucose-v8b.RDS")
phy_in <- phy.cpp

dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-chuckran-soil-glucose.rds")
sel <- which(phy.fxn@sam_data$Sample_name == "C0D3")
phy.fxn@sam_data$Run[sel] # "SRR9032615"
sel <- which(dat.cpd.collate$sample == "SRR9032615")
dat.cpd.collate <- dat.cpd.collate[-sel, ]
saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-excludeC0D3-chuckran-soil-glucose.rds")
dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-excludeC0D3-chuckran-soil-glucose.rds")

sum(sample_sums(phy_in)) # all: 867.1734 ; exclude C0D3: 800.5567
sample_sums(phy_in)
# exclude C0D3:
# SRR9032199 SRR9032202 SRR9032258 SRR9032259 SRR9032267 SRR9032300 SRR9032509 SRR9032510 SRR9032617 SRR9032694 SRR9032715 SRR9032716 
# 66.83239   66.76902   66.78190   66.70470   66.62648   66.74872   66.60375   66.57557   66.60223   66.70927   66.65317   66.94948 
# all:
# SRR9032199 SRR9032202 SRR9032258 SRR9032259 SRR9032267 SRR9032300 SRR9032509 SRR9032510 SRR9032615 SRR9032617 SRR9032694 SRR9032715 SRR9032716 
# 66.83239   66.76902   66.78190   66.70470   66.62648   66.74872   66.60375   66.57557   66.61675   66.60223   66.70927   66.65317   66.94948 

summary( sample_sums(phy_in) )
# exclude C0D3:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 66.58   66.62   66.71   66.71   66.77   66.95 

# all:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 66.58   66.62   66.70   66.71   66.77   66.95 

sd( sample_sums(phy_in) )
# exclude C0D3: 0.1101967
# all: 0.108834

max(taxa_sums(phy_in)) # # exclude C0D3: 95.68033 ; all: 103.6496


# don't rarefy - already in form of relative abundance %


## ordination plot
## PCoA + Bray-Curtis


## Already normalised to relative abundance

# # rarefy #1
# seed <- 1234
# r1 <- rarefy_even_depth(phy_in, sample.size = min(sample_sums(phy_in)),
#                         rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# min(taxa_sums(r1)) # 1
# sample_sums(r1) # all 3073
# ntaxa(r1) # 1014

r1 <- phy_in


### ORDINATION PLOT # # # # # # # # # # # # # # # 
### PCoA + Bray-Curtis

set.seed(1234)
ord <- ordinate(r1, "PCoA", "bray")

ord

str(r1@sam_data)

names(r1@sam_data)
# [1] "IMG_id"                    "Sample_name"               "Time_h"                    "JGI_analysis_project_type" "NCBI_BioProject_accession" "NCBI_BioSample_accession" 
# [7] "SRA_accession"             "No_reads"                  "Assembled_genome_size_bp"  "No_of_genes"               "No_of_scaffolds"           "N50_bp"                   
# [13] "GC_content"                "COG_database"              "Pfam_database"             "KEGG_database"             "Run"                       "fxn_sum_counts"           
# [19] "group_label"        


#saveRDS(r1, file = "r1-cpp-phyloseq-object-chuckran-soil-glucose-v8b.RDS")
saveRDS(r1, file = "r1-cpp-phyloseq-object-excludeC0D3-chuckran-soil-glucose-v8b.RDS")


p <- plot_ordination(r1, ord, type="samples", color="group_label")
#p <- plot_ordination(r1, ord, type="samples", color="Treatment_no_description", shape = "Treatment_no_description")
p <- plot_ordination(r1, ord, type="samples", color="group_label", label = "Sample_name")

## remove sample "C0D3" due to outlying beta diversity ordination position

p

str(p$data)

# x_lab <- p$labels$x
# y_lab <- p$labels$y

x_lab <- gsub(pattern = "Axis.1", replacement = "PCo1" , x = p$labels$x)
y_lab <- gsub(pattern = "Axis.2", replacement = "PCo2" , x =  p$labels$y)


names(p$data)
# [1] "Axis.1"                    "Axis.2"                    "IMG_id"                    "Sample_name"               "Time_h"                   
# [6] "JGI_analysis_project_type" "NCBI_BioProject_accession" "NCBI_BioSample_accession"  "SRA_accession"             "No_reads"                 
# [11] "Assembled_genome_size_bp"  "No_of_genes"               "No_of_scaffolds"           "N50_bp"                    "GC_content"               
# [16] "COG_database"              "Pfam_database"             "KEGG_database"             "Run"                       "fxn_sum_counts"           
# [21] "group_label"     


# cols <- c("0h" = "#8da0cb", 
#           "8h" = "#66c2a5", 
#           "24h" = "#fc8d62",
#           "48h" = "#e78ac3"
#           )

cols <- c("0h" = "#a6611a", 
          "8h" = "#dfc27d", 
          "24h" = "#80cdc1",
          "48h" = "#018571"
)

## add concave hulls to each group ...
# getting the convex hull of each unique point set
#https://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2
data_for_hulls <- p$data[ ,c("Sample_name","Axis.1","Axis.2","group_label")]
str(data_for_hulls)
levels(data_for_hulls$group_label)
# "0h"  "8h"  "24h" "48h"
find_hull <- function(data_for_hulls) data_for_hulls[chull(data_for_hulls$Axis.1, data_for_hulls$Axis.2), ]
hulls <- ddply(data_for_hulls, "group_label", find_hull)


label = "PERMANOVA: ~Time\nR^2 = 0.75, P = 0.001\nBeta-dispersion: P = 0.68"


pp <- ggplot(data=p$data, aes(x=Axis.1, y=Axis.2)) + # , colour=Sample_type__row_type x=NMDS1, y=NMDS2
  theme_bw() + 
  
  #geom_point(aes(colour=abbrev), size = 2) + # , alpha = 0.6
  #scale_color_manual(values=cols.bact, name = "Bacteria")+
  geom_point(aes(color=group_label), size = 2 ) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  scale_fill_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  geom_polygon(data=hulls,aes(group = group_label, fill = group_label), alpha = 0.175) + 
  
  #geom_text_repel(aes(label = abbrev), size = 3)+
  geom_text_npc(npcx = "left", npcy = "bottom", label = label, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE, 
  
  xlab(x_lab) + ylab(y_lab)+
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.title = element_text(size = rel(1.1)),
    
    legend.background = element_rect(fill = "transparent"),
    #legend.margin = margin(t = 0,r = 0,b = 0,l = 1,unit = "pt"),
    legend.margin = margin(t = 0,r = 0,b = 0,l = -5,unit = "pt"),
    legend.key.size = unit(0.8,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8))
  )

pp

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Beta-diversity-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 10, units = "cm", res = 600, compression = "lzw",type="cairo")



## PERMANOVA

# Calculate bray curtis distance matrix
set.seed(123)
bray <- phyloseq::distance(r1, method = "bray")
sampledf <- data.frame(sample_data(r1))
str(sampledf)

names(r1@sam_data)
# [1] "IMG_id"                    "Sample_name"               "Time_h"                    "JGI_analysis_project_type" "NCBI_BioProject_accession"
# [6] "NCBI_BioSample_accession"  "SRA_accession"             "No_reads"                  "Assembled_genome_size_bp"  "No_of_genes"              
# [11] "No_of_scaffolds"           "N50_bp"                    "GC_content"                "COG_database"              "Pfam_database"            
# [16] "KEGG_database"             "Run"                       "fxn_sum_counts"            "group_label"    

# Adonis test
set.seed(123)
adonis2(bray ~ group_label , data = sampledf)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = bray ~ group_label, data = sampledf)
#             Df   SumOfSqs      R2      F Pr(>F)    
# group_label  3 0.00031272 0.75094 8.0401  0.001 ***
# Residual     8 0.00010372 0.24906                  
# Total       11 0.00041644 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


beta <- betadisper(bray, sampledf$group_label)
set.seed(123)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
#           Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     3 7.1140e-06 2.3712e-06 0.5704    999  0.676
# Residuals  8 3.3254e-05 4.1568e-06  






### Alpha diversity

# Shannon Diversity Index
a_div <- plot_richness(r1, measures=c("Shannon")) #, "Simpson")) # Observed = Richness, but requires count data
a_div
# Shannon index emphasises richness, while Simpson index emphasises evenness

names(a_div$data)
# [1] "IMG_id"                    "Sample_name"               "Time_h"                    "JGI_analysis_project_type" "NCBI_BioProject_accession"
# [6] "NCBI_BioSample_accession"  "SRA_accession"             "No_reads"                  "Assembled_genome_size_bp"  "No_of_genes"              
# [11] "No_of_scaffolds"           "N50_bp"                    "GC_content"                "COG_database"              "Pfam_database"            
# [16] "KEGG_database"             "Run"                       "fxn_sum_counts"            "group_label"               "samples"                  
# [21] "variable"                  "value"                     "se"              


head(a_div$data)


# include Kruskal-Wallis test results - per below: 
#ktresult <- paste0("~Kruskal-Wallis:","~P == 7.6","~x ","~10^{-10}") # "\n",
ktresult <- "Kruskal-Wallis: P = 0.024"

#set.seed(123)
p <- ggplot(data=a_div$data, aes(x=group_label, y=value)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 ) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "Shannon diversity") +
  
  geom_text_npc(npcx = "middle", npcy = "bottom", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(5.525, 5.530, 5.537, 5.533), label = c("a","ab","b","ab"), size = 3.75)+
  
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Alpha-diversity-Shannon-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


names(a_div$data)
unique(a_div$data$variable) # "Shannon"
ktdat <- filter(a_div$data[ ,c("variable","value","group_label")], variable == "Shannon")
# Kruskal-Wallis test
kt <- kruskal.test( value ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  value by group_label
# Kruskal-Wallis chi-squared = 9.4615, df = 3, p-value = 0.02374
## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
unique( ktdat$group_label )
# [1] 8h  0h  24h 48h
# Levels: 0h < 8h < 24h < 48h
#pt <- dunnTest( value ~ group_label, data = ktdat, method = "bh") # Error in if (tmp$Eclass != "factor") { : the condition has length > 1
pt <- dunnTest( ktdat$value, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "h"   "24h" "48h" "8h" 
sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1     h      a         a 
# 4    8h     ab         ab
# 2   24h      b          b
# 3   48h     ab         ab
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "h"   "8h"  "24h" "48h"




####
#### Assess CPP in selected compounds
####

dat <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-excludeC0D3-chuckran-soil-glucose.rds")


## c) Glucose
sel.cpd <- which(df.comp$name == "D-Glucose")
this_var <- "Glucose"

df.comp[sel.cpd, ]
#             id    abbrev      name    form 
# 27    cpd00027     glc-D D-Glucose C6H12O6 
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6   

sel <- which(dat$cpd_id == "cpd00027")

head(dat[sel, ])

temp_dat <- dat[sel, ]

# include Kruskal-Wallis test results - per below: 
ktresult <- "Kruskal-Wallis: P = 0.022"


p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 ) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "Glucose CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "bottom", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(0.188, 0.1945, 0.2065, 0.2045), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-c-Glucose-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 9.6667, df = 3, p-value = 0.02162
## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
unique( ktdat$group_label )
# [1] 8h  0h  24h 48h
# Levels: 0h < 8h < 24h < 48h
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "h"   "24h" "48h" "8h" 
sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1     h      a         a 
# 4    8h     ab         ab
# 2   24h      b          b
# 3   48h     ab         ab
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "h"   "8h"  "24h" "48h"



## d) Cellulose
sel.cpd <- which(df.comp$name == "Cellulose")
this_var <- "Cellulose"

df.comp[sel.cpd, ]
#             id    abbrev      name      form  OC_ratio HC_ratio NC_ratio 
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2 0.8333333 1.666667        0

sel <- which(dat$cpd_id == "cpd11746")

head(dat[sel, ])
temp_dat <- dat[sel, ]

# # include Kruskal-Wallis test results - per below: 
# ktresult <- "Kruskal-Wallis: P = 0.022"

# include Kendall tau test results - per below:
ktresult <- "Kendall tau = 0.60, P = 0.011"

p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  #geom_boxplot(outlier.shape = NA)+
  geom_smooth(aes(group=1), method = "loess", alpha = 0.2)+
  
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "Cellulose CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  # annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(0.188, 0.1945, 0.2065, 0.2045), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 12.5,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-d-Cellulose-Glucose-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")

ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 6.2821, df = 3, p-value = 0.09867
# ## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
# unique( ktdat$group_label )
# # [1] 8h  0h  24h 48h
# # Levels: 0h < 8h < 24h < 48h
# pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
# pt
# pt$dtres
# pt <- pt$res
# pt
# cldList(comparison = pt$Comparison,
#         p.value    = pt$P.adj,
#         threshold  = 0.05)
# sig <- cldList(comparison = pt$Comparison,
#                p.value    = pt$P.adj,
#                threshold  = 0.05)
# str(sig)
# unique(sig$Group) # "h"   "24h" "48h" "8h" 
# sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
#                      ordered=TRUE )
# sig[ order(sig$Group), ]
# # Group Letter MonoLetter
# # 1     h      a         a 
# # 4    8h     ab         ab
# # 2   24h      b          b
# # 3   48h     ab         ab
# sig <- sig[ order(sig$Group), ]
# levels(sig$Group) # "h"   "8h"  "24h" "48h"


# Kendall tau correlation
df <- data.frame(x = as.numeric(temp_dat$group_label), y = temp_dat$cpd_rel_abun)
ktcor<- cor.test(x = df$x, y = df$y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  df$x and df$y
# z = 2.5584, p-value = 0.01052
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.6030227 





## e) CO2 - "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2

sel <- which(dat$cpd_id == "cpd00011")


head(dat[sel, ])
temp_dat <- dat[sel, ]

# include Kruskal-Wallis test results - per below: 
ktresult <- "Kruskal-Wallis: P = 0.04"

## include Kendall tau test results - per below:
#ktresult <- "Kendall tau = 0.60, P = 0.04"


p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "CO2 CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  #annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(0.188, 0.1945, 0.2065, 0.2045), label = c("a","ab","b","ab"), size = 3.75)+
  annotate("text", x = c(0.8, 1.8, 2.8, 3.8), y = c(0.918, 0.916, 0.907, 0.9145), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p


grid.text(label = "(e)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-e-CO2-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 8.2308, df = 3, p-value = 0.04148
## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
unique( ktdat$group_label )
# [1] 8h  0h  24h 48h
# Levels: 0h < 8h < 24h < 48h
pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "h"   "24h" "48h" "8h"
sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1     h      a         a 
# 4    8h     ab         ab
# 2   24h      b          b
# 3   48h     ab         ab
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "h"   "8h"  "24h" "48h"






## f) will be heatmap



## g) O2 - "Oxygen"
sel.cpd <- which(df.comp$name == "O2")
this_var <- "O2"

df.comp[sel.cpd, ]
#         id abbrev name form
# 7 cpd00007     o2   O2   O2

sel <- which(dat$cpd_id == "cpd00007")

head(dat[sel, ])
temp_dat <- dat[sel, ]

# # include Kruskal-Wallis test results - per below: 
# ktresult <- "Kruskal-Wallis: P = 0.022"

# include Kendall tau test results - per below:
ktresult <- "Kendall tau = 0.80, P = 0.0006"

p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  #geom_boxplot(outlier.shape = NA)+
  geom_smooth(aes(group=1), method = "loess", alpha = 0.2)+
  
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "O2 CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  # annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(0.188, 0.1945, 0.2065, 0.2045), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(g)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-g-O2-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp_dat
# Kruskal-Wallis test
kt <- kruskal.test( cpd_rel_abun ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  cpd_rel_abun by group_label
# Kruskal-Wallis chi-squared = 9.2564, df = 3, p-value = 0.02607

# ## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
# unique( ktdat$group_label )
# # [1] 8h  0h  24h 48h
# # Levels: 0h < 8h < 24h < 48h
# pt <- dunnTest( ktdat$cpd_rel_abun, ktdat$group_label, method = "bh")
# pt
# pt$dtres
# pt <- pt$res
# pt
# cldList(comparison = pt$Comparison,
#         p.value    = pt$P.adj,
#         threshold  = 0.05)
# sig <- cldList(comparison = pt$Comparison,
#                p.value    = pt$P.adj,
#                threshold  = 0.05)
# str(sig)
# unique(sig$Group) # "h"   "24h" "48h" "8h" 
# sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
#                      ordered=TRUE )
# sig[ order(sig$Group), ]
# # Group Letter MonoLetter
# # 1     h      a         a 
# # 4    8h     ab         ab
# # 2   24h      b          b
# # 3   48h     ab         ab
# sig <- sig[ order(sig$Group), ]
# levels(sig$Group) # "h"   "8h"  "24h" "48h"


# Kendall tau correlation
df <- data.frame(x = as.numeric(temp_dat$group_label), y = temp_dat$cpd_rel_abun)
ktcor<- cor.test(x = df$x, y = df$y, method = "kendall")
ktcor
# Kendall's rank correlation tau
# data:  df$x and df$y
# z = 3.4112, p-value = 0.0006467
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#       tau 
# 0.8040303 




# ## instead compare H2O ??
# 
# ## g) H2O - "Water"
# sel.cpd <- which(df.comp$name == "H2O")
# this_var <- "O2"
# 
# df.comp[sel.cpd, ]
# #         id abbrev name form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count
# # 1 cpd00001    h2o  H2O  H2O 
# 
# sel <- which(dat$cpd_id == "cpd00001")
# 
# head(dat[sel, ])
# temp_dat <- dat[sel, ]
# 
# 
# p <- ggplot(data=temp_dat, aes(x=group_label, y=cpd_rel_abun)) +
#   #theme_bw()+
#   theme_classic()+
#   geom_boxplot(outlier.shape = NA)+
#   #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
#   geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
#   scale_color_manual(values=cols.ecosystem_type, name = "Soil\nsource\nfor mice\ncages")+
#   
#   labs(x = NULL, y = "H2O CPP(%)") +
#   
#   theme(
#     #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
#     #legend.position = "none",
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.background = element_rect(fill = "transparent"),
#     legend.key.size = unit(0.9,"line"),
#     legend.title = element_text(size = rel(0.9)),
#     legend.text = element_text(size = rel(0.8)) ,
#     axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
#   )
# 
# p
# 
# 
# grid.text(label = "(g)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
# dev.print(tiff, file = paste0(workdir,"/plots/","CPP-g-H2O-Liu_2020_mice_soil_PRJNA542998-Hostremoval-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")




## h) 
# adenylate energy charge (AEC) indicates the energetic status of soil microorganisms
# the energy status of soilmicroorganisms was evaluated by determining AEC defined as: 
# AEC = (ATP + 0.5 × ADP) / (ATP + ADP + AMP)

sel.cpd <- which(df.comp$name == "ATP")
df.comp[sel.cpd, ]
#.        id abbrev name          form OC_ratio HC_ratio NC_ratio 
# 2 cpd00002    atp  ATP C10H13N5O13P3      1.3      1.3      0.5  
sel <- which(dat$cpd_id == "cpd00002")
head(dat[sel, ])
vals <- list()
vals[["ATP"]] <- dat[sel, ]
vals[["ATP"]]$sample

sel.cpd <- which(df.comp$name == "ADP")
df.comp[sel.cpd, ]
#         id abbrev name          form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio
# 8 cpd00008    adp  ADP C10H13N5O10P2        1      1.3      0.5 
sel <- which(dat$cpd_id == "cpd00008")
head(dat[sel, ])
vals[["ADP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["ADP"]]$sample ) # TRUE

sel.cpd <- which(df.comp$name == "AMP")
df.comp[sel.cpd, ]
#          id abbrev name        form OC_ratio HC_ratio NC_ratio
# 18 cpd00018    amp  AMP C10H12N5O7P      0.7      1.2      0.5 
sel <- which(dat$cpd_id == "cpd00018")
head(dat[sel, ])
vals[["AMP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["AMP"]]$sample ) # TRUE

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
AMP <- vals[["AMP"]]$cpd_rel_abun

AEC <- (ATP + 0.5*ADP) / (ATP + ADP + AMP)

temp <- cbind(dat[sel, ],data.frame(AEC=AEC))
head(temp)

# # include Kruskal-Wallis test results - per below: 
ktresult <- "Kruskal-Wallis: P = 0.02"

p <- ggplot(data=temp, aes(x=group_label, y=AEC)) +
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  labs(x = NULL, y = "AEC, from CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "bottom", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  annotate("text", x = c(1.2, 2.2, 3.2, 4.2), y = c(0.6686, 0.6694, 0.671, 0.6707), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(h)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-h-AEC-ratio-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp
# Kruskal-Wallis test
kt <- kruskal.test( AEC ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  AEC by group_label
# Kruskal-Wallis chi-squared = 9.6667, df = 3, p-value = 0.02162
## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
unique( ktdat$group_label )
# [1] 8h  0h  24h 48h
# Levels: 0h < 8h < 24h < 48h
pt <- dunnTest( ktdat$AEC, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "h"   "24h" "48h" "8h"
sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1     h      a         a 
# 4    8h     ab         ab
# 2   24h      b          b
# 3   48h     ab         ab
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "h"   "8h"  "24h" "48h"





## i) ATP / ADP

# use vectors for ATP and ADP from above

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
#AMP <- vals[["AMP"]]$cpd_rel_abun

ATP_ADP_ratio <- ATP/ADP

temp <- cbind(dat[sel, ],data.frame(ATP_ADP_ratio=ATP_ADP_ratio))

# # include Kruskal-Wallis test results - per below: 
ktresult <- "Kruskal-Wallis: P = 0.02"


p <- ggplot(data=temp, aes(x=group_label, y=ATP_ADP_ratio)) +
  theme_classic()+
  geom_boxplot(outlier.shape = NA)+
  #geom_point(aes(color=group_label), size = 2 , alpha = 0.6) + # , alpha = 0.6
  geom_jitter(aes(color=group_label), width =0.2, height = 0, size = 2 , alpha = 0.6) + # , alpha = 0.6
  scale_color_manual(values=cols, name = "Time\nafter\nglucose\naddition")+
  
  labs(x = NULL, y = "ATP/ADP, from CPP(%)") +
  
  geom_text_npc(npcx = "middle", npcy = "top", label = ktresult, size = 3.25 , lineheight = 0.85 )+ # parse = TRUE,
  annotate("text", x = c(0.8, 1.8, 2.8, 3.8), y = c(1.665, 1.663, 1.6485, 1.652), label = c("a","ab","b","ab"), size = 3.75)+
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    #legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(size = rel(1.1))
    #axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(i)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-i-ATP_ADP_ratio-chuckran-soil-glucose-v8b.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")


ktdat <- temp
# Kruskal-Wallis test
kt <- kruskal.test( ATP_ADP_ratio ~ group_label, data = ktdat ) # Kruskal Wallis test
kt
# Kruskal-Wallis rank sum test
# data:  ATP_ADP_ratio by group_label
# Kruskal-Wallis chi-squared = 9.6667, df = 3, p-value = 0.02162
## Dunn Test uses factor vector or non-numeric vector that can be coerced to a factor vector
unique( ktdat$group_label )
# [1] 8h  0h  24h 48h
# Levels: 0h < 8h < 24h < 48h
pt <- dunnTest( ktdat$ATP_ADP_ratio, ktdat$group_label, method = "bh")
pt
pt$dtres
pt <- pt$res
pt
cldList(comparison = pt$Comparison,
        p.value    = pt$P.adj,
        threshold  = 0.05)
sig <- cldList(comparison = pt$Comparison,
               p.value    = pt$P.adj,
               threshold  = 0.05)
str(sig)
unique(sig$Group) # "h"   "24h" "48h" "8h"
sig$Group <- factor( sig$Group, levels=c("h", "8h", "24h", "48h"  ),
                     ordered=TRUE )
sig[ order(sig$Group), ]
# Group Letter MonoLetter
# 1     h      a         a 
# 4    8h     ab         ab
# 2   24h      b          b
# 3   48h     ab         ab
sig <- sig[ order(sig$Group), ]
levels(sig$Group) # "h"   "8h"  "24h" "48h"




## f) HEATMAP

dim(dat) # 101832      6
head(dat)
#     cpd_id     sample cpd_rel_abun log10_abun group group_label
# 1 cpd25681 SRR9032199 0.0004249709 -3.3716409     8          8h
# 2 cpd02597 SRR9032199 0.0191206785 -1.7184967     8          8h
# 3 cpd24620 SRR9032199 0.0005762276 -3.2394059     8          8h
# 4 cpd00001 SRR9032199 5.0011883693  0.6990732     8          8h
# 5 cpd01501 SRR9032199 0.0131099817 -1.8823979     8          8h
# 6 cpd00851 SRR9032199 0.0108064371 -1.9663175     8          8h

# p<- ggplot(dat, aes(x = sample, y = cpd_id, fill = log10_abun)) + # ggplot(long_df, aes(x = column_name, y = row_id, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_minimal()
# p


#library(pheatmap)

dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ sample', value.var = "log10_abun")
#dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ abbrev', value.var = "cpd_rel_abun")
row.names(dat.wide) <- dat.wide[ ,1]
dat.wide <- dat.wide[ ,-1]

# delete compounds with minimal variation
#sel <- which(dat.wide == 0) # 

# calculate row standard deviation?
row_sd <- apply(dat.wide, 1, sd, na.rm = TRUE)
hist(row_sd); summary(row_sd)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.009939 0.021772 0.073687 0.052266 1.712593


dim(dat.wide) # 8486   12
length(row_sd) # 8486
#low_sd_rows <- df[row_sd < 0.5, ]

quantile(row_sd, probs = 0.5)
# cpd_rel_abun
# 50% 
# 0.0217721 
sel <- which(row_sd > quantile(row_sd, probs = 0.5) ) # 4232

dat.wide.select <- dat.wide[sel, ]

colnames(dat.wide.select)
identical( colnames(dat.wide.select) , rownames(phy_in@sam_data) ) # TRUE

# col.df <- data.frame(sample=colnames(dat.wide.select), group_label=phy_in@sam_data$group_label)
# row.names(col.df) <- col.df$sample

col.df <- data.frame(Time=phy_in@sam_data$group_label)
row.names(col.df) <- colnames(dat.wide.select)

cols
# 0h        8h       24h       48h 
# "#a6611a" "#dfc27d" "#80cdc1" "#018571" 

ann_colors = list(
  Time = c( `0h` = "#a6611a", `8h` = "#dfc27d", `24h` = "#80cdc1", `48h` = "#018571"  )
)

#pheatmap(as.matrix(dat.wide))
#pheatmap(as.matrix(dat.wide), scale = "row", show_rownames = FALSE)

#pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 10)

pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 10, show_colnames = FALSE, annotation_col = col.df , annotation_colors = ann_colors )


grid.text(label = "(f)", x = unit(0.033, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=18, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-f-HEATMAP-chuckran-soil-glucose-v8b.tiff"), width = 18, height = 22, units = "cm", res = 500, compression = "lzw",type="cairo")


#-------------------------


##########################
##########################
##########################
##########################

## Supplementary validation dataset

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures
#    prepare for data download from SRA
#-------------------------

# Reference WGS culture
# PRJNA622674_NIBSC_WGS_cultures
# Amos, G.C.A. et al. Developing standards for the microbiome field. Microbiome 8, 98 (2020).
# In this study, we describe the development of the first reference reagents produced by the National Institute for Biological Standards and Control (NIBSC) for microbiome analysis by next-generation sequencing.
# All sequencing data generated through this study is publicly available through the NCBI Sequence Read Archive upon publication (NCBI Bioproject ID PRJNA622674).


sradat <- read_excel(path = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/SraRunTable_PRJNA622674_NIBSC_Standards.xlsx", sheet = 1, range = "A1:AL56")
sradat <- as.data.frame(sradat)

unique(sradat$LibrarySource) # "METAGENOMIC" "GENOMIC"
sel <- which(sradat$LibrarySource == "GENOMIC")
sradat <- sradat[sel, ]

row.names(sradat) <- sradat$Run

str(sradat)
# 'data.frame':	20 obs. of  38 variables:
# $ Run                           : chr  "SRR11487909" "SRR11487910" "SRR11487911" "SRR11487912" ...
# $ Assay Type                    : chr  "WGS" "WGS" "WGS" "WGS" ...
# $ AvgSpotLen                    : num  302 298 298 298 298 298 302 302 298 298 ...
# $ Bases                         : num  1.28e+09 1.07e+09 1.29e+09 1.51e+09 1.31e+09 ...
# $ BioProject                    : chr  "PRJNA622674" "PRJNA622674" "PRJNA622674" "PRJNA622674" ...
# $ BioSample                     : chr  "SAMN14524783" "SAMN14524782" "SAMN14524781" "SAMN14524780" ...
# $ BioSampleModel                : chr  "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" ...
# $ Bytes                         : num  5.00e+08 4.00e+08 4.85e+08 5.74e+08 4.95e+08 ...
# $ Center Name                   : chr  "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" ...
# $ Collection_Date               : POSIXct, format: "2017-01-10" "2017-01-10" "2017-01-10" "2017-01-10" ...
# $ Consent                       : chr  "public" "public" "public" "public" ...
# $ DATASTORE filetype            : chr  "fastq,sra,run.zq" "fastq,run.zq,sra" "fastq,run.zq,sra" "sra,run.zq,fastq" ...
# $ DATASTORE provider            : chr  "s3,ncbi,gs" "ncbi,gs,s3" "gs,ncbi,s3" "s3,gs,ncbi" ...
# $ DATASTORE region              : chr  "gs.us-east1,s3.us-east-1,ncbi.public" "ncbi.public,s3.us-east-1,gs.us-east1" "ncbi.public,s3.us-east-1,gs.us-east1" "ncbi.public,s3.us-east-1,gs.us-east1" ...
# $ Experiment                    : chr  "SRX8063932" "SRX8063931" "SRX8063930" "SRX8063929" ...
# $ geo_loc_name_country          : chr  "United Kingdom" "United Kingdom" "United Kingdom" "United Kingdom" ...
# $ geo_loc_name_country_continent: chr  "Europe" "Europe" "Europe" "Europe" ...
# $ geo_loc_name                  : chr  "United Kingdom" "United Kingdom" "United Kingdom" "United Kingdom" ...
# $ Instrument                    : chr  "NextSeq 500" "NextSeq 500" "NextSeq 500" "NextSeq 500" ...
# $ isolation_source              : chr  "culture" "culture" "culture" "culture" ...
# $ lat_lon                       : chr  "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" ...
# $ Library Name                  : chr  "294_294-Sample8_S66" "306_Sample_16_S16" "306_Sample_22_S22" "316_7089_S36" ...
# $ LibraryLayout                 : chr  "PAIRED" "PAIRED" "PAIRED" "PAIRED" ...
# $ LibrarySelection              : chr  "RANDOM" "RANDOM" "RANDOM" "RANDOM" ...
# $ LibrarySource                 : chr  "GENOMIC" "GENOMIC" "GENOMIC" "GENOMIC" ...
# $ Organism                      : chr  "Ruminococcus gauvreauii" "Roseburia intestinalis" "Roseburia hominis" "Prevotella melaninogenica" ...
# $ Lifestyle                     : chr  "anaerobic acetate producer" "thrives on fibre-rich diet" "thrives on fibre-rich diet" "upper respiratory tract, opportunistic" ...
# $ Platform                      : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
# $ ReleaseDate                   : chr  "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" ...
# $ create_date                   : chr  "2020-04-07T04:11:00Z" "2020-04-07T04:32:00Z" "2020-04-07T04:40:00Z" "2020-04-07T04:18:00Z" ...
# $ version                       : num  1 1 1 1 1 1 1 1 1 1 ...
# $ Sample Name                   : chr  "DSM_19829" "DSM_14610" "DSM_16839" "DSM_7089" ...
# $ SRA Study                     : chr  "SRP255413" "SRP255413" "SRP255413" "SRP255413" ...
# $ HOST                          : chr  "DSM_19829" "DSM_14610" "DSM_16839" "DSM_7089" ...
# $ isolate                       : chr  "missing" "missing" "missing" "missing" ...
# $ sample_type                   : chr  "missing" "missing" "missing" "missing" ...
# $ strain                        : chr  "missing" "missing" "missing" "missing" ...
# $ sub_species                   : chr  NA NA NA NA ...

sradat$Organism
# [1] "Ruminococcus gauvreauii"                "Roseburia intestinalis"                 "Roseburia hominis"                      "Prevotella melaninogenica"             
# [5] "Segatella copri"                        "Parabacteroides distasonis"             "Lactobacillus gasseri"                  "Anaerobutyricum hallii"                
# [9] "Escherichia coli"                       "Collinsella aerofaciens"                "Clostridium butyricum"                  "Blautia wexlerae"                      
# [13] "Bifidobacterium longum subsp. longum"   "Bifidobacterium longum subsp. infantis" "Bacteroides uniformis"                  "Faecalibacterium prausnitzii"          
# [17] "Bacteroides thetaiotaomicron"           "Anaerostipes hadrus"                    "Alistipes finegoldii"                   "Akkermansia muciniphila"     

runlist <- read.table(file = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/PRJNA622674_NIBSC_WGS_cultures-run-list.txt", header = FALSE )
runlist <- as.character(runlist[ ,1])

identical(runlist, row.names(sradat)) # TRUE


sradat.select <- sradat

dim(sradat.select) # 20 38


saveRDS(object = sradat.select, file = "sradat.select.PRJNA622674_NIBSC_WGS_cultures.rds")


# file used for SRA download was: "PRJNA622674_NIBSC_WGS_cultures-run-list.txt"


#-------------------------

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - read in superfocus - fxn potential outputs
#-------------------------

#saveRDS(object = sradat.select, file = "sradat.select.PRJNA622674_NIBSC_WGS_cultures.rds")

sradat.select <- readRDS("sradat.select.PRJNA622674_NIBSC_WGS_cultures.rds")

sampid <- sradat.select$Run # 20

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487909"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487910"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487911"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487912"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487913"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487915"


names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR11487909 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487909" 
# SRR11487910 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487910" 
# SRR11487911 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487911" 
# SRR11487912 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487912" 
# SRR11487913 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487913" 
# SRR11487915 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/Validation-Reference-Samples/3_fxn_superfocus_copy/superfocus_out_SRR11487915"  


length(results_dirs) # 20

sel <- which(names(results_dirs) %in% sampid) # qty 20
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 20

# check identical order
identical(sampid, names(results_dirs)) # TRUE
#identical(sort(sampid), sort(names(results_dirs))) #
length(results_dirs) # 20
length(sampid) # 20

# # reset sampid to remove missing sample
# sampid <- names(results_dirs)
# identical(sampid, names(results_dirs)) # TRUE


# In this data one Run corresponds to a single Sample_ID !!!

# collate results into a long-format table

sfx.long <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,percent_abun=NA)

for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  # [1] "Subsystem.Level.1"                                                                         
  # [2] "Subsystem.Level.2"                                                                         
  # [3] "Subsystem.Level.3"                                                                         
  # [4] "Function"                                                                                  
  # [5] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq"  
  # [6] "X.scratch.pawsey1216.cliddicoat.ft2d_chn.2b_clean_hostremoval.SRR341581_non_host.1.fastq.."
  
  
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # [5] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [6] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [7] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq"
  # [8] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq"
  # [9] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [10] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H2THFBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  # [11] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L001_R1.good.fastq.." # this is %
  # [12] "X.scratch.user.lidd0026.ami_2_fastp_qc.12465_1_PE_550bp_BASE_UNSW_H3WYJBCXX_TAATGCGC.TAATCTTA_L002_R1.good.fastq.." # this is %
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.fastq..$", x = names(tab)) # for single (unpaired) reads
  if (length(sel.col.percent)>1) {
    tab$percent_abun <- apply(X = tab[ ,sel.col.percent], MARGIN = 1, FUN = mean )
  } else {
    tab$percent_abun <- tab[ ,sel.col.percent]
  }
  
  # sum(tab$percent_abun) # 100
  # mean(tab$percent_abun) # 0.004338583
  
  names(sfx.long) # "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"    "percent_abun"
  # names(tab)
  # [1] "Subsystem.Level.1"
  # [2] "Subsystem.Level.2"
  # [3] "Subsystem.Level.3"
  # [4] "Function"
  # ...
  # [13] "sampid"
  # [14] "percent_abun"
  
  tab <- tab[ ,c("sampid","Subsystem.Level.1","Subsystem.Level.2","Subsystem.Level.3","Function","percent_abun")]
  names(tab) <- names(sfx.long)
  
  sfx.long <- rbind(sfx.long,tab)
  
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}


head(sfx.long)
# remove empty 1st row
sfx.long <- sfx.long[-1, ]
dim(sfx.long) # 58177     6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                           subsys_L3
# 2 SRR11487909 Amino Acids and Derivatives                            -                 Amino acid racemase
# 3 SRR11487909 Amino Acids and Derivatives                            -                 Amino acid racemase
# 4 SRR11487909 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 5 SRR11487909 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 6 SRR11487909 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 7 SRR11487909 Amino Acids and Derivatives Alanine, serine, and glycine                Alanine biosynthesis
# fxn percent_abun
# 2                           UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 2.351992e-05
# 3 UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 2.351992e-05
# 4                                                                                            Creatinine_amidohydrolase_(EC_3.5.2.10) 8.072037e-02
# 5                                                                                           Hydantoinase/oxoprolinase_family_protein 1.223036e-03
# 6                                                                               N-methylhydantoinase_(ATP-hydrolyzing)_(EC_3.5.2.14) 6.115179e-04
# 7                                                                                                      Alanine_racemase_(EC_5.1.1.1) 1.030408e-01

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 12537    21

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 12537
length(unique(full_fxn_names)) # 12537

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Aspartate_racemase_(EC_5.1.1.13)" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Diaminopimelate_epimerase_(EC_5.1.1.7)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Glutamate_racemase_(EC_5.1.1.3)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR11487909"  "SRR11487910"  "SRR11487911"  "SRR11487912"  "SRR11487913"  "SRR11487915"  "SRR11487916"  "SRR11487917"  "SRR11487918" 
# [11] "SRR11487919"  "SRR11487920"  "SRR11487921"  "SRR11487922"  "SRR11487923"  "SRR11487924"  "SRR11487926"  "SRR11487927"  "SRR11487928"  "SRR11487929" 
# [21] "SRR11487930" 

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR11487909" "SRR11487910" "SRR11487911" "SRR11487912" "SRR11487913" "SRR11487915"

length(sampid) # 20

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # TRUE
#identical(sort(names(sfx.wide)), sort(as.character(sampid))) #

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 12537     4

length(unique(tax.fxn$subsys_L1)) # 34
length(unique(tax.fxn$subsys_L2)) # 170
length(unique(tax.fxn$subsys_L3)) # 932
length(unique(tax.fxn$fxn)) # 6839


# # # #

## gather Function count data??
sfx.long.count <- data.frame(sampleID=NA, subsys_L1=NA, subsys_L2=NA, subsys_L3=NA,fxn=NA,count_abun=NA)
length(sampid) # 20
for (i in 1:length(sampid)) {
  #i<-1
  this_samp <- sampid[i]
  sel.folder <- grep(pattern = this_samp, x = results_dirs)
  this_folder <- results_dirs[sel.folder]
  #tab1 <- read_excel(path = paste0(this_folder,"/output_all_levels_and_function.xlsx"), skip = 4, col_names = TRUE)
  tab <- read.csv(file = paste0(this_folder,"/output_all_levels_and_function.xls"), sep = "\t", skip = 4 )
  # names(tab)
  tab$sampid <- this_samp
  names(tab)
  tab <- tab[,c(7,1,2,3,4,5)] # this time capture 'count' data
  names(tab) <- names(sfx.long.count)
  sfx.long.count <- rbind(sfx.long.count,tab)
  print(paste0("completed ",i," - sample ID: ",sampid[i]))
}
head(sfx.long.count)
# remove empty 1st row
sfx.long.count <- sfx.long.count[-1, ]
sum(sfx.long.count$count_abun) # 28619375 = 28,619,375
sfx.long.count$full_fxn_tax <- paste0(sfx.long.count$subsys_L1,"___", sfx.long.count$subsys_L2,"___", sfx.long.count$subsys_L3,"___", sfx.long.count$fxn)
head(sfx.long.count)
sfx.wide.count <- dcast(sfx.long.count, formula = full_fxn_tax ~ sampleID, value.var = "count_abun")
dim(sfx.wide.count) # 12537    21
sel.na <- which(is.na(sfx.wide.count),arr.ind = TRUE)
sfx.wide.count[sel.na] <- 0
sum(colSums(sfx.wide.count[,-1])) # 28619375
hist(colSums(sfx.wide.count[,-1]))
mean(colSums(sfx.wide.count[,-1])) # 1430969
sd(colSums(sfx.wide.count[,-1])) # 374306.3

summary(colSums(sfx.wide.count[,-1]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 946423 1062596 1446382 1430969 1622033 2379931 
length(unique(sfx.long.count$subsys_L1)) # 34

fxn_sum_counts <- colSums(sfx.wide.count[,-1])

# # # #


#-------------------------

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - get into Phyloseq object
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 12537     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 12537    20

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12537 taxa and 20 samples ]
# tax_table()   Taxonomy Table:    [ 12537 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR11487909" "SRR11487910" "SRR11487911" "SRR11487912" "SRR11487913" "SRR11487915" "SRR11487916" "SRR11487917" "SRR11487918" "SRR11487919" "SRR11487920"
# [12] "SRR11487921" "SRR11487922" "SRR11487923" "SRR11487924" "SRR11487926" "SRR11487927" "SRR11487928" "SRR11487929" "SRR11487930"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

head(row.names(sradat.select))

samp <- sradat.select

dim(samp) # 20 38

head(row.names(samp)) # 

#row.names(samp) <- samp$Run

identical(row.names(samp), sample_names(phy)) # TRUE

#length(row.names(samp)) # 102
#length(sample_names(phy)) # 101
#sel <- which(row.names(samp) %in% sample_names(phy)) # 101
#samp2 <- samp[sample_names(phy), ]

identical(row.names(samp), names(fxn_sum_counts)) # TRUE

samp$fxn_sum_counts <- fxn_sum_counts


SAMP <- sample_data(samp)



### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12537 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 39 sample variables ]
# tax_table()   Taxonomy Table:    [ 12537 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Aspartate_racemase_(EC_5.1.1.13)"             
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Diaminopimelate_epimerase_(EC_5.1.1.7)"       
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Glutamate_racemase_(EC_5.1.1.3)"         


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"

saveRDS(object = phy, file = "phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS")


head(phy@sam_data)

# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 20
NonZeroFxns
# SRR11487909 SRR11487910 SRR11487911 SRR11487912 SRR11487913 SRR11487915 SRR11487916 SRR11487917 SRR11487918 SRR11487919 SRR11487920 SRR11487921 SRR11487922 SRR11487923 
# 3880        2994        2933        2070        2205        2880        2170        2627        7208        2097        3536        2873        2423        2548 
# SRR11487924 SRR11487926 SRR11487927 SRR11487928 SRR11487929 SRR11487930 
# 3050        2389        3346        2579        2174        2195 

mean(NonZeroFxns) # 2908.85
sd(NonZeroFxns) # 1131.336


#-------------------------

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - COPY of R code to run CPP steps on HPC
#    1) build reaction search - get reactions and compounds
#    2) get cpd rel abun per sample
#    3) collate compounds for each sample
#-------------------------

# # # # # # # # # # # # #
# #
# # R script for cpp3d
# # - build reaction search in parallel - get_reactions & compounds
# # - get cpd rel abun per sample
# # - collate_compounds
# #
# # For study - PRJNA622674_NIBSC_WGS_cultures
# # Craig Liddicoat - Flinders University
# # Running on Pawsey Setonix
# # # # # # # # # # # # #
# 
# # Add a new path
# .libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
#             "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))
# 
# R.Version()
# 
# # load packages
# #library(readxl); packageVersion("readxl")
# library(parallel); packageVersion("parallel")
# library(doParallel); packageVersion("doParallel")
# library(dplyr); packageVersion("dplyr")
# library(stringr); packageVersion("stringr")
# library(phyloseq); packageVersion("phyloseq") # '1.44.0'
# 
# message("\n# establish folders and input files")
# 
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis'")
# workdir <- "/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis/working"
# 
# message("\nthis_study <- '-nibsc-ref-cultures-pawsey'")
# this_study <- "-nibsc-ref-cultures-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS")
# 
# 
# subsys.lut <- readRDS("subsys.lut.RDS")
# rxns.lut <- readRDS("rxns.lut.RDS")
# rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
# compounds.lut <- readRDS("compounds.lut.RDS")
# 
# 
# 
# message("\n### 1) build reaction search in parallel - get_reactions & compounds")
# message("\n# # # # # # # # # #")
# message("\ndf.tax <- as.data.frame(phy@tax_table)")
# df.tax <- as.data.frame(phy@tax_table)
# message("\nhead(row.names(df.tax))")
# head(row.names(df.tax))
# message("\ndim(df.tax)")
# dim(df.tax)
# 
# 
# get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
#   
#   rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#   rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#   
#   sub1 <- df.tax$subsys_L1[i]
#   sub2 <- df.tax$subsys_L2[i]
#   sub3 <- df.tax$subsys_L3[i]
#   
#   fxn.temp <- df.tax$fxn[i]
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   
#   # store results corresponding to each Superfocus row
#   fxn.list <- list()
#   fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
#   
#   # check for multiple functions/reactions?
#   flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
#   flag2 <- grepl(pattern = "_@_", x = fxn.temp)
#   if (!any(flag1,flag2)==TRUE) {
#     # no multiples
#     fxns <- fxn.temp
#   } else if (flag1==TRUE) {
#     fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
#   } else {
#     fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
#   }
#   # remove underscores
#   ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
#   
#   # process each fxn & store attributes
#   df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
#   
#   # Identify '/' separators with no '_'  ??
#   
#   for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
#     #f<-1
#     #f<-2
#     f.in <- fxns[f]
#     
#     # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
#     full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
#     full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
#     
#     ## data cleaning
#     
#     # trim off '_#' and '_##' tags
#     trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
#     if (!is.na(trim_nchar) & length(trim_nchar)==1) {
#       f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
#     }
#     
#     # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
#     f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
#     
#     #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
#     #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
#     
#     full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
#     full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
#     
#     sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
#     
#     ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
#     if (length(sel.rx)>=1) {
#       df.fxns$matching_method[f] <- "Exact hierachy match"
#       df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#       
#     } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
#       # search by EC id if present
#       
#       f.in <- fxns[f] # this goes back to string with brackets for EC
#       ## LOOK FOR MULTIPLE ECs ?
#       
#       how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
#       
#       ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
#       #class(ECs)
#       ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
#       ECs.collapse <- paste0(ECs, collapse = "|")
#       
#       sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
#       
#       if (length(how_many_ECs)==0 | length(ECs)==0) {
#         # there was a glitch, database typo, or some error in identifying the EC number
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       } else if (length(sel.rx)>=1) {
#         # combined EC hits identified
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
#         # treat EC hits individually
#         sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
#         # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
#         sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
#         
#         df.fxns$matching_method[f] <- "EC number"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else {
#         # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       # END EC matching
#       
#       
#     } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
#       ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
#       # 1. 'reactions' table by name: rxns.lut$name
#       # 2. 'reactions' table by aliases: rxns.lut$aliases
#       # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
#       # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
#       
#       if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
#         # 1a - exact match - rxns.lut$name
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
#         #rxns.lut$name[sel.rx]
#         df.fxns$matching_method[f] <- "Matched Reactions name"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
#         # 2a - exact match - rxns.lut$aliases
#         sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
#         #rxns.lut$aliases[sel.rx]
#         #rxns.lut$name[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Reactions aliases"
#         df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
#         # 3a - exact match - subsys.lut$Role
#         sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
#         #subsys.lut$Role[sel.rx]
#         #subsys.lut$Reaction[sel.rx]
#         
#         df.fxns$matching_method[f] <- "Matched Subsytem role"
#         df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
#         
#       } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
#         # 4a - exact match - rxn_pathways.lut$External_rxn_name
#         sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
#         
#         df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
#         df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
#         
#         
#       } else {
#         df.fxns$matching_method[f] <- "No match found"
#         df.fxns$rxns[f] <- NA
#         
#       }
#       
#       ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
#       
#       
#     } # END function - reaction search
#     
#     #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
#     
#     #print(paste0("completed fxn ", f))
#     
#     
#     ## now investigate these reactions ...
#     # Reactions lookup table: 
#     # - "equation": Definition of reaction expressed using compound IDs and after protonation
#     # Compounds lookup table:
#     # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
#     #df.fxns
#     
#     
#     #if (df.fxns$matching_method == "No match found") {
#     if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
#       
#       df.Rxns <- NA
#       df.Compounds <- NA
#       
#     } else { # reaction(s) were identified
#       
#       # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
#       f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
#       
#       df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
#                             rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
#       
#       for (r in 1:dim(df.Rxns)[1]) {
#         #r<-1
#         #this_rxn <- "rxn00004"
#         this_rxn <- df.Rxns$rxn_id[r]
#         sel <- which(rxns.lut$id == this_rxn)
#         ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
#         ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
#         ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
#         
#         # extract compound info
#         
#         #df.Rxns$rxn_eqn[r]
#         #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
#         #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
#         
#         ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
#         # 5 23 43 61
#         # 6 25 46 65 83
#         
#         ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
#         # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
#         
#         if (length(compds)>=1) {
#           
#           df.Rxns$compds[r] <- paste0(compds, collapse = ";")
#           
#           ## get compound coefficients?
#           start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
#           end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
#           ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
#           
#           df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
#           
#           # get formulas of compounds
#           
#           formx <-filter(compounds.lut, id %in% compds )
#           row.names(formx) <- formx$id
#           ( formx.char <- formx[compds, ]$formula )
#           # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
#           # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
#           # "C7H7O7" "H2O"    "C7H5O6"
#           df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
#           
#           ( compd.names <- formx[compds, ]$name )
#           # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
#           
#           temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
#                                           cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
#                                           #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
#           )
#           
#         } else {
#           # No specified reaction equation or chemical formula info
#           df.Rxns$compds[r] <- NA
#           df.Rxns$compd_coef[r] <- NA
#           df.Rxns$chem_formx[r] <- NA
#           
#           temp.df.Compounds <- NA
#           
#         }
#         
#         if (r==1) { df.Compounds <- temp.df.Compounds }
#         
#         if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
#         
#         # clean up - if there are additional reactions?
#         temp.df.Compounds <- NA
#         
#       } # END loop for r - rxn_id's per f/f.in
#       
#     } # END else loop when reactions identified
#     
#     # store results corresponding to each sub-reaction of each Superfocus row
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
#     
#     if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
#     fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
#     
#     
#   } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
#   
#   saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
#   
# } # END function to be run in parallel for each superfocus row
# 
# 
# # # # # # # # # # # # # # # # # # #
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system (may need to use other alternative for Windows)
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
#   get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# message("\n(num_results_files <- dim(df.tax)[1])")
# (num_results_files <- dim(df.tax)[1])
# 
# # assemble all compound data outputs
# # start with blank row
# 
# df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )
# 
# for (i in 1:num_results_files) {
#   fxn.superfocus.rowlabel <- row.names(df.tax)[i]
#   temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
#   
#   f_no <- length( temp[[1]][["compounds"]] )
#   
#   for (f in 1:f_no) {
#     #f<-2
#     # only add non-NA results
#     if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
#       
#       df.temp <- temp[[1]][["compounds"]][[f]]
#       ok <- complete.cases(df.temp)
#       df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
#       df.out <- rbind(df.out,df.temp)
#     }
#   }
#   print(paste0("added df ",i," of ",num_results_files ))
#   
# }
# 
# 
# message("\nstr(df.out)")
# str(df.out)
# 
# 
# saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))
# 
# # remove NA first row
# message("\nhead(df.out)")
# head(df.out)
# 
# df.out <- df.out[-1, ]
# 
# message("\ndim(df.out)")
# dim(df.out)
# 
# 
# message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")
# 
# df.out$cpd_molar_prop_norm <- NA
# 
# message("\nlength(unique(df.out$superfocus_fxn))")
# length(unique(df.out$superfocus_fxn))
# 
# message("\nphy")
# phy
# 
# message("\n% of functions represented - with compound information")
# 100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))
# 
# 
# fxns_found <- unique(df.out$superfocus_fxn)
# 
# for (k in 1:length(fxns_found)) {
#   #k<-1
#   this_fxn <- fxns_found[k]
#   sel <- which(df.out$superfocus_fxn == this_fxn)
#   
#   sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
#   # calculate 
#   
#   df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
#   
#   print(paste0("completed ",k))
#   
# }
# 
# message("\nsum(df.out$cpd_molar_prop_norm)")
# sum(df.out$cpd_molar_prop_norm)
# 
# message("\nsample_sums(phy)")
# sample_sums(phy)
# 
# message("\ngetwd()")
# getwd()
# 
# saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))
# 
# 
# 
# message("\n### 2) get cpd rel abun per sample")
# message("\n# # # # # # # # # #")
# 
# 
# df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
# message("\ndim(df.OTU)")
# dim(df.OTU)
# 
# 
# get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
#   
#   this_samp <- sample_names(phy_in)[i]
#   df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
#   
#   dat.cpd$sample <- this_samp
#   
#   dat.cpd$cpd_rel_abun_norm <- NA
#   
#   fxns_all <- row.names(df.OTU)
#   
#   for (k in 1:length(fxns_all)) {
#     #k<-1
#     this_fxn <- fxns_all[k]
#     sel <- which(dat.cpd$superfocus_fxn == this_fxn)
#     
#     if (length(sel)>=1) {
#       dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
#       
#     }
#   } # END rel abun values for all relevant functions added
#   
#   saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   
# } # END
# 
# 
# no_forks <- 8
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
# head(dat)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
#   dat <- rbind(dat, temp)
#   print(paste0("completed ",i))
# }
# 
# 
# saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )
# 
# rm(temp)
# 
# message("\nstr(dat)")
# str(dat)
# 
# message("\nsum(dat$cpd_rel_abun_norm)")
# sum(dat$cpd_rel_abun_norm)
# 
# message("\naverage functional relative abundance per sample")
# message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# 
# message("\nnames(dat)")
# names(dat)
# 
# message("\nlength(unique(dat$cpd_id))")
# length(unique(dat$cpd_id))
# 
# 
# 
# 
# message("\n### 3) collate_compounds within each sample")
# message("\n# # # # # # # # # #")
# 
# 
# unique_cpd <- unique(dat$cpd_id)
# samp_names <- sample_names(phy)
# 
# 
# collate_compounds <- function(dat.cpd, unique_cpd, samp) {
#   #i<-1
#   #samp = samp_names[i]
#   #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
#   
#   this_samp <- samp
#   
#   cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
#   
#   for (c in 1:length(unique_cpd)) {
#     #c<-1
#     this_cpd <- unique_cpd[c]
#     sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
#     
#     if (length(sel.cpd) >=1) {
#       cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
#     }
#     
#   } # END all compounds
#   
#   saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
# } # END
# 
# 
# 
# no_forks <- 4
# 
# # this makes clusters on Unix-like system
# cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
# registerDoParallel(cl)
# 
# foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
#   collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])
# 
# stopCluster(cl)
# 
# 
# message("\n## assemble results")
# 
# # output 1
# i<-1
# this_samp <- sample_names(phy)[i]
# dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
# head(dat.cpd.collate)
# 
# for ( i in 2:length(sample_names(phy)) ) {
#   this_samp <- sample_names(phy)[i]
#   temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
#   
#   dat.cpd.collate <- rbind(dat.cpd.collate, temp)
#   
#   print(paste0("completed ",i))
# }
# 
# 
# message("\nstr(dat.cpd.collate)")
# str(dat.cpd.collate)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)")
# sum(dat.cpd.collate$cpd_rel_abun)
# 
# message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# 
# saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))
# 
# # END


#-------------------------

# PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - COPY of OUTOUTS from R code after running CPP steps on HPC
#-------------------------

# $platform
# [1] "x86_64-pc-linux-gnu"
# 
# $arch
# [1] "x86_64"
# 
# $os
# [1] "linux-gnu"
# 
# $system
# [1] "x86_64, linux-gnu"
# 
# $status
# [1] ""
# 
# $major
# [1] "4"
# 
# $minor
# [1] "4.1"
# 
# $year
# [1] "2024"
# 
# $month
# [1] "06"
# 
# $day
# [1] "14"
# 
# $`svn rev`
# [1] "86737"
# 
# $language
# [1] "R"
# 
# $version.string
# [1] "R version 4.4.1 (2024-06-14)"
# 
# $nickname
# [1] "Race for Your Life"
# 
# [1] ‘4.4.1’
# Loading required package: foreach
# Loading required package: iterators
# [1] ‘1.0.17’
# 
# Attaching package: ‘dplyr’
# 
# The following objects are masked from ‘package:stats’:
#   
#   filter, lag
# 
# The following objects are masked from ‘package:base’:
#   
#   intersect, setdiff, setequal, union
# 
# [1] ‘1.1.4’
# [1] ‘1.5.2’
# [1] ‘1.46.0’
# 
# # establish folders and input files
# 
# workdir <- '/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis/working'
# 
# this_study <- '-nibsc-ref-cultures-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS')
# 
# ### 1) build reaction search in parallel - get_reactions & compounds
# 
# # # # # # # # # # #
# 
# df.tax <- as.data.frame(phy@tax_table)
# 
# head(row.names(df.tax))
# [1] "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"
# 
# dim(df.tax)
# [1] 12537     4
# [[1]]
# NULL
# 
# ...
# 
# 
# [[12536]]
# NULL
# 
# [[12537]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 12537
# [1] "added df 1 of 12537"
# [1] "added df 2 of 12537"
# [1] "added df 3 of 12537"
# ...
# 
# 
# [1] "added df 12536 of 12537"
# [1] "added df 12537 of 12537"
# 
# str(df.out)
# 'data.frame':	316702 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_1" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn00283" "rxn00283" "rxn19085" ...
# $ cpd_id        : chr  NA "cpd00035" "cpd00117" "cpd00035" ...
# $ cpd_name      : chr  NA "L-Alanine" "D-Alanine" "L-Alanine" ...
# $ cpd_form      : chr  NA "C3H7NO2" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                         f__in   rxn_id   cpd_id  cpd_name
# 1           <NA> NA                          <NA>     <NA>     <NA>      <NA>
#   2          fxn_1  1 Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035 L-Alanine
# 3          fxn_1  1 Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117 D-Alanine
# 4          fxn_1  1 Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035 L-Alanine
# 5          fxn_1  1 Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117 D-Alanine
# 6          fxn_1  1 Alanine racemase (EC 5.1.1.1) rxn38030 cpd00035 L-Alanine
# cpd_form cpd_molar_prop
# 1     <NA>             NA
# 2  C3H7NO2              1
# 3  C3H7NO2              1
# 4  C3H7NO2              1
# 5  C3H7NO2              1
# 6  C3H7NO2              1
# 
# dim(df.out)
# [1] 316701      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 6806
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12537 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 39 sample variables ]
# tax_table()   Taxonomy Table:    [ 12537 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.28731
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 6805"
# [1] "completed 6806"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 6806
# 
# sample_sums(phy)
# SRR11487909 SRR11487910 SRR11487911 SRR11487912 SRR11487913 SRR11487915 
# 100         100         100         100         100         100 
# SRR11487916 SRR11487917 SRR11487918 SRR11487919 SRR11487920 SRR11487921 
# 100         100         100         100         100         100 
# SRR11487922 SRR11487923 SRR11487924 SRR11487926 SRR11487927 SRR11487928 
# 100         100         100         100         100         100 
# SRR11487929 SRR11487930 
# 100         100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/PRJNA622674_NIBSC_WGS_cultures/cpp_analysis"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 12537    20
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# ...
# 
# 
# 
# [[20]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                         f__in   rxn_id   cpd_id  cpd_name
# 2          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035 L-Alanine
# 3          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117 D-Alanine
# 4          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035 L-Alanine
# 5          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117 D-Alanine
# 6          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn38030 cpd00035 L-Alanine
# 7          fxn_1 1 Alanine racemase (EC 5.1.1.1) rxn38030 cpd00117 D-Alanine
# cpd_form cpd_molar_prop cpd_molar_prop_norm      sample cpd_rel_abun_norm
# 2  C3H7NO2              1           0.1666667 SRR11487909                 0
# 3  C3H7NO2              1           0.1666667 SRR11487909                 0
# 4  C3H7NO2              1           0.1666667 SRR11487909                 0
# 5  C3H7NO2              1           0.1666667 SRR11487909                 0
# 6  C3H7NO2              1           0.1666667 SRR11487909                 0
# 7  C3H7NO2              1           0.1666667 SRR11487909                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# [1] "completed 18"
# [1] "completed 19"
# [1] "completed 20"
# 
# str(dat)
# 'data.frame':	6334020 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_1" "fxn_1" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn00283" "rxn00283" "rxn19085" "rxn19085" ...
# $ cpd_id             : chr  "cpd00035" "cpd00117" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "L-Alanine" "D-Alanine" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C3H7NO2" "C3H7NO2" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.167 0.167 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR11487909" "SRR11487909" "SRR11487909" "SRR11487909" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 1274.486
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 63.72432
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6192
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# [[1]]
# NULL
# ...
# 
# [[20]]
# NULL
# 
# 
# ## assemble results
# cpd_id      sample cpd_rel_abun
# 1 cpd00035 SRR11487909 2.238630e-01
# 2 cpd00117 SRR11487909 2.140593e-01
# 3 cpd00041 SRR11487909 2.150706e-01
# 4 cpd00320 SRR11487909 1.941327e-05
# 5 cpd00504 SRR11487909 5.496867e-02
# 6 cpd00516 SRR11487909 9.451592e-02
# [1] "completed 2"
# [1] "completed 3"
# ...
# [1] "completed 19"
# [1] "completed 20"
# 
# str(dat.cpd.collate)
# 'data.frame':	123840 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00035" "cpd00117" "cpd00041" "cpd00320" ...
# $ sample      : chr  "SRR11487909" "SRR11487909" "SRR11487909" "SRR11487909" ...
# $ cpd_rel_abun: num  2.24e-01 2.14e-01 2.15e-01 1.94e-05 5.50e-02 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 1274.486
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 63.72432
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries

#-------------------------

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - continue CPP analysis
#-------------------------

"NIBSC_WGS_cultures_PRJNA622674"

phy.fxn <- readRDS('phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS')
phy.fxn
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12537 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 39 sample variables ]
# tax_table()   Taxonomy Table:    [ 12537 taxa by 4 taxonomic ranks ]

phy <- phy.fxn

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/NIBSC_WGS_cultures_PRJNA622674/cpp_analysis/dat.cpd.collate-all-samps-cpp3d--nibsc-ref-cultures-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	123840 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00035" "cpd00117" "cpd00041" "cpd00320" ...
# $ sample      : chr  "SRR11487909" "SRR11487909" "SRR11487909" "SRR11487909" ...
# $ cpd_rel_abun: num  2.24e-01 2.14e-01 2.15e-01 1.94e-05 5.50e-02 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000055 0.010291 0.001266 7.799974 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -4.2589    -Inf -2.8974  0.8921 


# log10 abun
dat.cpd.collate$log10_abun <- dat.cpd.collate$cpd_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value
subsel.zero <- which(dat.cpd.collate$log10_abun == 0) # qty 52948
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(dat.cpd.collate$log10_abun[ -subsel.zero ])
  dat.cpd.collate$log10_abun[ subsel.zero ] <- zero_replace
}
dat.cpd.collate$log10_abun <- log10(dat.cpd.collate$log10_abun)

hist(dat.cpd.collate$log10_abun); summary( dat.cpd.collate$log10_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.3074 -8.3074 -4.2589 -5.3181 -2.8974  0.8921 


# make group variable from sample name

dat.cpd.collate$Organism <- NA
dat.cpd.collate$abbrev <- NA

samp <- as(phy.fxn@sam_data, "data.frame")
names(samp)
# [1] "Run"                            "Assay.Type"                     "AvgSpotLen"                     "Bases"                         
# [5] "BioProject"                     "BioSample"                      "BioSampleModel"                 "Bytes"                         
# [9] "Center.Name"                    "Collection_Date"                "Consent"                        "DATASTORE.filetype"            
# [13] "DATASTORE.provider"             "DATASTORE.region"               "Experiment"                     "geo_loc_name_country"          
# [17] "geo_loc_name_country_continent" "geo_loc_name"                   "Instrument"                     "isolation_source"              
# [21] "lat_lon"                        "Library.Name"                   "LibraryLayout"                  "LibrarySelection"              
# [25] "LibrarySource"                  "Organism"                       "Lifestyle"                      "Platform"                      
# [29] "ReleaseDate"                    "create_date"                    "version"                        "Sample.Name"                   
# [33] "SRA.Study"                      "HOST"                           "isolate"                        "sample_type"                   
# [37] "strain"                         "sub_species"                    "fxn_sum_counts"     

samp$Organism
# [1] "Ruminococcus gauvreauii"                "Roseburia intestinalis"                 "Roseburia hominis"                     
# [4] "Prevotella melaninogenica"              "Segatella copri"                        "Parabacteroides distasonis"            
# [7] "Lactobacillus gasseri"                  "Anaerobutyricum hallii"                 "Escherichia coli"                      
# [10] "Collinsella aerofaciens"                "Clostridium butyricum"                  "Blautia wexlerae"                      
# [13] "Bifidobacterium longum subsp. longum"   "Bifidobacterium longum subsp. infantis" "Bacteroides uniformis"                 
# [16] "Faecalibacterium prausnitzii"           "Bacteroides thetaiotaomicron"           "Anaerostipes hadrus"                   
# [19] "Alistipes finegoldii"                   "Akkermansia muciniphila"  

samp$abbrev <- c(
  "R. gauvreauii"      ,   "R. intestinalis"        ,   "R. hominis"     ,         
  "P. melaninogenica"  ,   "S. copri"               ,   "P. distasonis"  ,          
  "L. gasseri"         ,   "A. hallii"              ,   "E. coli"        ,          
  "C. aerofaciens"     ,   "C. butyricum"           ,   "B. wexlerae"    ,          
  "B. longum ssp longum",   "B. longum ssp infantis", "B. uniformis"     ,          
  "F. prausnitzii"     ,   "B. thetaiotaomicron"    ,   "A. hadrus"      ,          
  "A. finegoldii"      ,   "A. muciniphila"  
)

identical( phy@sam_data$Run , samp$Run ) # TRUE
identical( sample_names(phy), samp$Run ) # TRUE


#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$Organism[sel] <- as.character( samp$Organism[i] )
  dat.cpd.collate$abbrev[sel] <- as.character( samp$abbrev[i] )
  print(paste0("completed ", i))
}


head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-NIBSC_WGS_cultures_PRJNA622674.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-NIBSC_WGS_cultures_PRJNA622674.rds")


str(dat.cpd.collate)
# 'data.frame':	123840 obs. of  6 variables:
# $ cpd_id      : chr  "cpd00035" "cpd00117" "cpd00041" "cpd00320" ...
# $ sample      : chr  "SRR11487909" "SRR11487909" "SRR11487909" "SRR11487909" ...
# $ cpd_rel_abun: num  2.24e-01 2.14e-01 2.15e-01 1.94e-05 5.50e-02 ...
# $ log10_abun  : num  -0.65 -0.669 -0.667 -4.712 -1.26 ...
# $ Organism    : chr  "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" ...
# $ abbrev      : chr  "R. gauvreauii" "R. gauvreauii" "R. gauvreauii" "R. gauvreauii" ...


length( unique(dat.cpd.collate$cpd_id) ) # 6192
6192*20 # 123840


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id      sample cpd_rel_abun log10_abun                Organism        abbrev
# 1 cpd00035 SRR11487909 2.238630e-01 -0.6500177 Ruminococcus gauvreauii R. gauvreauii
# 2 cpd00117 SRR11487909 2.140593e-01 -0.6694659 Ruminococcus gauvreauii R. gauvreauii
# 3 cpd00041 SRR11487909 2.150706e-01 -0.6674190 Ruminococcus gauvreauii R. gauvreauii
# 4 cpd00320 SRR11487909 1.941327e-05 -4.7119014 Ruminococcus gauvreauii R. gauvreauii
# 5 cpd00504 SRR11487909 5.496867e-02 -1.2598848 Ruminococcus gauvreauii R. gauvreauii
# 6 cpd00516 SRR11487909 9.451592e-02 -1.0244950 Ruminococcus gauvreauii R. gauvreauii

dim(data_in) # 123840      6

unique_samps <- unique(data_in$sample)

no_compounds <- numeric(length = length(unique_samps))
sample_sum_relabun <- numeric(length = length(unique_samps))

for (i in 1:length(unique_samps)) {
  #i<-1
  this_samp <- unique_samps[i]
  sel <- which(data_in$sample == this_samp)
  
  values <- data_in$cpd_rel_abun[sel]
  values <- values[values > 0]
  
  no_compounds[i] <- length( values )
  sample_sum_relabun[i] <- sum(values)
  print(paste0("completed ",i))
}

mean(no_compounds) # 3544.6
sd(no_compounds) #  546.0895

mean(sample_sum_relabun) # 63.72432
sd(sample_sum_relabun) # 5.23906

length(unique(data_in$cpd_id)) # 6192

#-------------------------

#### PRJNA622674_NIBSC_WGS_cultures - nibsc-ref-cultures - continue CPP analysis
#    CPP - get into phyloseq object
#    beta diversity
#    alpha diversity
#    response in selected compounds: glucose, cellulose, CO2, O2, AEC, ATP/ADP
#    heatmap of scaled CPP
#-------------------------

phy.fxn <- readRDS('phy-phyloseq-fxn-PRJNA622674_NIBSC_WGS_cultures.RDS')
phy.fxn
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 12537 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 39 sample variables ]
# tax_table()   Taxonomy Table:    [ 12537 taxa by 4 taxonomic ranks ]

str(samp) # this is sample metadata from above
# 'data.frame':	20 obs. of  40 variables:
# $ Run                           : chr  "SRR11487909" "SRR11487910" "SRR11487911" "SRR11487912" ...
# $ Assay.Type                    : chr  "WGS" "WGS" "WGS" "WGS" ...
# $ AvgSpotLen                    : num  302 298 298 298 298 298 302 302 298 298 ...
# $ Bases                         : num  1.28e+09 1.07e+09 1.29e+09 1.51e+09 1.31e+09 ...
# $ BioProject                    : chr  "PRJNA622674" "PRJNA622674" "PRJNA622674" "PRJNA622674" ...
# $ BioSample                     : chr  "SAMN14524783" "SAMN14524782" "SAMN14524781" "SAMN14524780" ...
# $ BioSampleModel                : chr  "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" "Microbe\\, viral or environmental" ...
# $ Bytes                         : num  5.00e+08 4.00e+08 4.85e+08 5.74e+08 4.95e+08 ...
# $ Center.Name                   : chr  "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" "NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL" ...
# $ Collection_Date               : POSIXct, format: "2017-01-10" "2017-01-10" "2017-01-10" "2017-01-10" ...
# $ Consent                       : chr  "public" "public" "public" "public" ...
# $ DATASTORE.filetype            : chr  "fastq,sra,run.zq" "fastq,run.zq,sra" "fastq,run.zq,sra" "sra,run.zq,fastq" ...
# $ DATASTORE.provider            : chr  "s3,ncbi,gs" "ncbi,gs,s3" "gs,ncbi,s3" "s3,gs,ncbi" ...
# $ DATASTORE.region              : chr  "gs.us-east1,s3.us-east-1,ncbi.public" "ncbi.public,s3.us-east-1,gs.us-east1" "ncbi.public,s3.us-east-1,gs.us-east1" "ncbi.public,s3.us-east-1,gs.us-east1" ...
# $ Experiment                    : chr  "SRX8063932" "SRX8063931" "SRX8063930" "SRX8063929" ...
# $ geo_loc_name_country          : chr  "United Kingdom" "United Kingdom" "United Kingdom" "United Kingdom" ...
# $ geo_loc_name_country_continent: chr  "Europe" "Europe" "Europe" "Europe" ...
# $ geo_loc_name                  : chr  "United Kingdom" "United Kingdom" "United Kingdom" "United Kingdom" ...
# $ Instrument                    : chr  "NextSeq 500" "NextSeq 500" "NextSeq 500" "NextSeq 500" ...
# $ isolation_source              : chr  "culture" "culture" "culture" "culture" ...
# $ lat_lon                       : chr  "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" "51.6884 N 0.2409 W" ...
# $ Library.Name                  : chr  "294_294-Sample8_S66" "306_Sample_16_S16" "306_Sample_22_S22" "316_7089_S36" ...
# $ LibraryLayout                 : chr  "PAIRED" "PAIRED" "PAIRED" "PAIRED" ...
# $ LibrarySelection              : chr  "RANDOM" "RANDOM" "RANDOM" "RANDOM" ...
# $ LibrarySource                 : chr  "GENOMIC" "GENOMIC" "GENOMIC" "GENOMIC" ...
# $ Organism                      : chr  "Ruminococcus gauvreauii" "Roseburia intestinalis" "Roseburia hominis" "Prevotella melaninogenica" ...
# $ Lifestyle                     : chr  "anaerobic acetate producer" "thrives on fibre-rich diet" "thrives on fibre-rich diet" "upper respiratory tract, opportunistic" ...
# $ Platform                      : chr  "ILLUMINA" "ILLUMINA" "ILLUMINA" "ILLUMINA" ...
# $ ReleaseDate                   : chr  "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" "2020-04-30T00:00:00Z" ...
# $ create_date                   : chr  "2020-04-07T04:11:00Z" "2020-04-07T04:32:00Z" "2020-04-07T04:40:00Z" "2020-04-07T04:18:00Z" ...
# $ version                       : num  1 1 1 1 1 1 1 1 1 1 ...
# $ Sample.Name                   : chr  "DSM_19829" "DSM_14610" "DSM_16839" "DSM_7089" ...
# $ SRA.Study                     : chr  "SRP255413" "SRP255413" "SRP255413" "SRP255413" ...
# $ HOST                          : chr  "DSM_19829" "DSM_14610" "DSM_16839" "DSM_7089" ...
# $ isolate                       : chr  "missing" "missing" "missing" "missing" ...
# $ sample_type                   : chr  "missing" "missing" "missing" "missing" ...
# $ strain                        : chr  "missing" "missing" "missing" "missing" ...
# $ sub_species                   : chr  NA NA NA NA ...
# $ fxn_sum_counts                : num  981165 1027843 1238611 1510358 1306958 ...
# $ abbrev                        : chr  "R. gauvreauii" "R. intestinalis" "R. hominis" "P. melaninogenica" ...
temp <- samp


dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-NIBSC_WGS_cultures_PRJNA622674.rds")

str(dat.cpd.collate)
# 'data.frame':	123840 obs. of  6 variables:
# $ cpd_id      : chr  "cpd00035" "cpd00117" "cpd00041" "cpd00320" ...
# $ sample      : chr  "SRR11487909" "SRR11487909" "SRR11487909" "SRR11487909" ...
# $ cpd_rel_abun: num  2.24e-01 2.14e-01 2.15e-01 1.94e-05 5.50e-02 ...
# $ log10_abun  : num  -0.65 -0.669 -0.667 -4.712 -1.26 ...
# $ Organism    : chr  "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" "Ruminococcus gauvreauii" ...
# $ abbrev      : chr  "R. gauvreauii" "R. gauvreauii" "R. gauvreauii" "R. gauvreauii" ...

data_in <- dat.cpd.collate
str(data_in)

length( unique(data_in$cpd_id) ) # 6192
length( unique(data_in$cpd_id[data_in$cpd_rel_abun > 0]) ) # 6192
length( unique(data_in$sample) ) # 20



### get data into phyloseq object ...

head(data_in)
#     cpd_id      sample cpd_rel_abun log10_abun                Organism        abbrev
# 1 cpd00035 SRR11487909 2.238630e-01 -0.6500177 Ruminococcus gauvreauii R. gauvreauii
# 2 cpd00117 SRR11487909 2.140593e-01 -0.6694659 Ruminococcus gauvreauii R. gauvreauii
# 3 cpd00041 SRR11487909 2.150706e-01 -0.6674190 Ruminococcus gauvreauii R. gauvreauii
# 4 cpd00320 SRR11487909 1.941327e-05 -4.7119014 Ruminococcus gauvreauii R. gauvreauii
# 5 cpd00504 SRR11487909 5.496867e-02 -1.2598848 Ruminococcus gauvreauii R. gauvreauii
# 6 cpd00516 SRR11487909 9.451592e-02 -1.0244950 Ruminococcus gauvreauii R. gauvreauii


df.wide <- dcast(data_in, formula = sample + Organism + abbrev ~ cpd_id , value.var = "cpd_rel_abun" )

df.wide[1:5, 1:10]
#        sample                  Organism            abbrev cpd00001 cpd00002  cpd00003  cpd00004  cpd00005  cpd00006   cpd00007
# 1 SRR11487909   Ruminococcus gauvreauii     R. gauvreauii 4.051303 3.261371 0.9195612 0.8903872 0.5486324 0.5534446 0.06650808
# 2 SRR11487910    Roseburia intestinalis   R. intestinalis 5.375296 2.645240 0.3944702 0.3699013 0.3642720 0.3660460 0.05465324
# 3 SRR11487911         Roseburia hominis        R. hominis 4.291290 2.684270 0.4495062 0.4162309 0.4337927 0.4375048 0.06443941
# 4 SRR11487912 Prevotella melaninogenica P. melaninogenica 4.456092 2.468731 0.5158637 0.4852961 0.4067525 0.4079656 0.06529177
# 5 SRR11487913           Segatella copri          S. copri 5.227585 2.382694 0.5669993 0.5349707 0.5286283 0.5304102 0.06039110


# save group variable
#samp <- df.wide[ ,1:2]
samp <- df.wide[ ,1:3]
row.names(samp) <- samp$sample

# transpose
#df.wide <- t(df.wide[ ,-2]) # minus 'group' column
df.wide <- t(df.wide[ ,-c(2,3)]) # minus 'Organism' and 'abbrev' columns

head(df.wide)

samp_names <- df.wide[1, ]
tax_names <- row.names(df.wide[-1, ])
head(tax_names) # "cpd00001" "cpd00002" "cpd00003" "cpd00004" "cpd00005" "cpd00006"
otu.df <- df.wide[-1, ] # remove sample labels in 1st row
# this is necessary to create numeric matrix

colnames(otu.df) <- samp_names

# convert OTU table to matrix
class(otu.df) # "matrix" "array"
#otu.df <- as.matrix(otu.df)

# convert to numeric matrix
# https://stackoverflow.com/questions/20791877/convert-character-matrix-into-numeric-matrix
otu.df <- apply(otu.df, 2, as.numeric)

rownames(otu.df) # NULL
dim(otu.df) #  6192   20
rownames(otu.df) <- tax_names

## Create 'otuTable'
#  otu_table - Works on any numeric matrix.
#  You must also specify if the species are rows or columns
OTU <- otu_table(otu.df, taxa_are_rows = TRUE)


# # convert Taxonomy table to matrix

tax <- data.frame(cpd_id = tax_names)
row.names(tax) <- tax_names

tax <- as.matrix(tax)

identical( row.names(otu.df), row.names(tax) ) # TRUE


## Create 'taxonomyTable'
#  tax_table - Works on any character matrix.
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
TAX <- tax_table(tax)


## Create a phyloseq object, merging OTU & TAX tables
phy.cpp = phyloseq(OTU, TAX)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6192 taxa and 20 samples ]
# tax_table()   Taxonomy Table:    [ 6192 taxa by 1 taxonomic ranks ]


sample_names(phy.cpp)
# [1] "SRR11487909" "SRR11487910" "SRR11487911" "SRR11487912" "SRR11487913" "SRR11487915" "SRR11487916" "SRR11487917" "SRR11487918" "SRR11487919" "SRR11487920"
# [12] "SRR11487921" "SRR11487922" "SRR11487923" "SRR11487924" "SRR11487926" "SRR11487927" "SRR11487928" "SRR11487929" "SRR11487930"

#identical(sample_names(phy.cpp), samp$sample) # TRUE

#identical(sample_names(phy.cpp), sradat.select2$Run) # TRUE
identical(sample_names(phy.cpp), temp$Run) # TRUE

#row.names(sradat.select2) <- sradat.select2$Run
identical( row.names(temp), sample_names(phy.cpp) ) # TRUE

#samp <- sradat.select2
samp <- temp

# row.names need to match sample_names() from phyloseq object
#row.names(samp) <- samp$sample
#identical(row.names(samp), samp$sample) # TRUE



### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

SAMP <- sample_data(samp)


### Combine SAMPDATA into phyloseq object
phy.cpp <- merge_phyloseq(phy.cpp, SAMP)
phy.cpp
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6192 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 40 sample variables ]
# tax_table()   Taxonomy Table:    [ 6192 taxa by 1 taxonomic ranks ]

phy.cpp@sam_data

min(taxa_sums(phy.cpp)) #  2.028369e-07

saveRDS(object = phy.cpp, file = "phy.cpp-cleaned-NIBSC_WGS_cultures_PRJNA622674-v8.RDS")
phy.cpp <- readRDS("phy.cpp-cleaned-NIBSC_WGS_cultures_PRJNA622674-v8.RDS")

phy_in <- phy.cpp

sum(sample_sums(phy_in)) # 1274.486
sample_sums(phy_in)
# SRR11487909 SRR11487910 SRR11487911 SRR11487912 SRR11487913 SRR11487915 SRR11487916 SRR11487917 SRR11487918 SRR11487919 SRR11487920 SRR11487921 SRR11487922 
# 66.32605    63.91754    63.79122    60.29156    65.68034    58.15840    61.36657    67.74977    52.48048    71.64224    60.97792    65.52249    70.92441 
# SRR11487923 SRR11487924 SRR11487926 SRR11487927 SRR11487928 SRR11487929 SRR11487930 
# 69.30488    59.65057    64.82713    53.20990    68.01563    62.92003    67.72927 

summary( sample_sums(phy_in) )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.48   60.81   64.37   63.72   67.73   71.64 

sd( sample_sums(phy_in) )
# 5.23906

max(taxa_sums(phy_in)) # 133.4995


# don't rarefy - already in form of relative abundance %


## ordination plot
## PCoA + Bray-Curtis


## Already normalised to relative abundance

# # rarefy #1
# seed <- 1234
# r1 <- rarefy_even_depth(phy_in, sample.size = min(sample_sums(phy_in)),
#                         rngseed = seed, replace = FALSE, trimOTUs = TRUE, verbose = TRUE)
# min(taxa_sums(r1)) # 1
# sample_sums(r1) # all 3073
# ntaxa(r1) # 1014

r1 <- phy_in


### ORDINATION PLOT # # # # # # # # # # # # # # # 
### PCoA + Bray-Curtis

set.seed(1234)
ord <- ordinate(r1, "PCoA", "bray")

ord

str(r1@sam_data)

names(r1@sam_data)
#  [1] "Run"                            "Assay.Type"                     "AvgSpotLen"                     "Bases"                         
# [5] "BioProject"                     "BioSample"                      "BioSampleModel"                 "Bytes"                         
# [9] "Center.Name"                    "Collection_Date"                "Consent"                        "DATASTORE.filetype"            
# [13] "DATASTORE.provider"             "DATASTORE.region"               "Experiment"                     "geo_loc_name_country"          
# [17] "geo_loc_name_country_continent" "geo_loc_name"                   "Instrument"                     "isolation_source"              
# [21] "lat_lon"                        "Library.Name"                   "LibraryLayout"                  "LibrarySelection"              
# [25] "LibrarySource"                  "Organism"                       "Lifestyle"                      "Platform"                      
# [29] "ReleaseDate"                    "create_date"                    "version"                        "Sample.Name"                   
# [33] "SRA.Study"                      "HOST"                           "isolate"                        "sample_type"                   
# [37] "strain"                         "sub_species"                    "fxn_sum_counts"                 "abbrev"    



#shapes <- c(0:18,25)

n <- 20
#set.seed(123)
#palette <- distinctColorPalette(n)
#pie(rep(1, n), col=palette)


#library(RColorBrewer); packageVersion("RColorBrewer") # ‘1.1.3’

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(123)
pie(rep(1,n), col=sample(col_vector, n))

set.seed(123)
palette <- sample(col_vector, n)

cols.bact <-c(
  palette
)

names(cols.bact) <- samp$abbrev

cols.bact
# R. gauvreauii        R. intestinalis             R. hominis      P. melaninogenica               S. copri          P. distasonis 
# "#CCEBC5"              "#FFFF33"              "#E6AB02"              "#80B1D3"              "#E6F5C9"              "#FF7F00" 
# L. gasseri              A. hallii                E. coli         C. aerofaciens           C. butyricum            B. wexlerae 
# "#FFF2AE"              "#BC80BD"              "#CAB2D6"              "#8DA0CB"              "#B15928"              "#BEBADA" 
# B. longum ssp longum B. longum ssp infantis           B. uniformis         F. prausnitzii    B. thetaiotaomicron              A. hadrus 
# "#1B9E77"              "#FBB4AE"              "#E5D8BD"              "#666666"              "#6A3D9A"              "#BF5B17" 
# A. finegoldii         A. muciniphila 
# "#FCCDE5"              "#B3B3B3" 



saveRDS(r1, file = "r1-cpp-phyloseq-object-NIBSC_WGS_cultures_PRJNA622674-v8.RDS")


p <- plot_ordination(r1, ord, type="samples", color="abbrev")
#p <- plot_ordination(r1, ord, type="samples", color="Treatment_no_description", shape = "Treatment_no_description")
p

str(p$data)

# x_lab <- p$labels$x
# y_lab <- p$labels$y

x_lab <- gsub(pattern = "Axis.1", replacement = "PCo1" , x = p$labels$x)
y_lab <- gsub(pattern = "Axis.2", replacement = "PCo2" , x =  p$labels$y)


names(p$data)
# [1] "Axis.1"                         "Axis.2"                         "Run"                            "Assay.Type"                    
# [5] "AvgSpotLen"                     "Bases"                          "BioProject"                     "BioSample"                     
# [9] "BioSampleModel"                 "Bytes"                          "Center.Name"                    "Collection_Date"               
# [13] "Consent"                        "DATASTORE.filetype"             "DATASTORE.provider"             "DATASTORE.region"              
# [17] "Experiment"                     "geo_loc_name_country"           "geo_loc_name_country_continent" "geo_loc_name"                  
# [21] "Instrument"                     "isolation_source"               "lat_lon"                        "Library.Name"                  
# [25] "LibraryLayout"                  "LibrarySelection"               "LibrarySource"                  "Organism"                      
# [29] "Lifestyle"                      "Platform"                       "ReleaseDate"                    "create_date"                   
# [33] "version"                        "Sample.Name"                    "SRA.Study"                      "HOST"                          
# [37] "isolate"                        "sample_type"                    "strain"                         "sub_species"                   
# [41] "fxn_sum_counts"                 "abbrev"   




pp <- ggplot(data=p$data, aes(x=Axis.1, y=Axis.2)) + # , colour=Sample_type__row_type x=NMDS1, y=NMDS2
  theme_bw() + 
  
  #geom_point(aes(colour=abbrev), size = 2) + # , alpha = 0.6
  #scale_color_manual(values=cols.bact, name = "Bacteria")+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  geom_text_repel(aes(label = abbrev), size = 3)+
  
  xlab(x_lab) + ylab(y_lab)+
  
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #plot.title = element_text(size = rel(1.1)),
    
    legend.background = element_rect(fill = "transparent"),
    #legend.margin = margin(t = 0,r = 0,b = 0,l = 1,unit = "pt"),
    legend.margin = margin(t = 0,r = 0,b = 0,l = -5,unit = "pt"),
    legend.key.size = unit(0.8,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8))
  )

pp

grid.text(label = "(a)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Beta-diversity-NIBSC_WGS_cultures_PRJNA622674.tiff"), width = 13, height = 10, units = "cm", res = 600, compression = "lzw",type="cairo")





### Alpha diversity

plot_richness(r1, measures=c("Observed")) 
# error in evaluating the argument 'x' in selecting a method for function 't': function accepts only integers (counts)

# Shannon Diversity Index
a_div <- plot_richness(r1, measures=c("Shannon")) #, "Simpson")) # Observed = Richness, but requires count data
a_div
# Shannon index emphasises richness, while Simpson index emphasises evenness

names(a_div$data)
#  "sample_id"      "Tube_name"      "Rep"            "Time"           "Carbon_source"  "Treatment"      "fxn_sum_counts" "samples"        "variable"       "value"          "se"  

head(a_div$data)
# Run Assay.Type AvgSpotLen      Bases  BioProject    BioSample                    BioSampleModel     Bytes
# 1 SRR11487909        WGS        302 1277976118 PRJNA622674 SAMN14524783 Microbe\\, viral or environmental 500016789
# 2 SRR11487910        WGS        298 1070105782 PRJNA622674 SAMN14524782 Microbe\\, viral or environmental 400303753
# 3 SRR11487911        WGS        298 1285275192 PRJNA622674 SAMN14524781 Microbe\\, viral or environmental 485084010
# 4 SRR11487912        WGS        298 1508737346 PRJNA622674 SAMN14524780 Microbe\\, viral or environmental 573672806
# 5 SRR11487913        WGS        298 1305936128 PRJNA622674 SAMN14524779 Microbe\\, viral or environmental 495163310
# 6 SRR11487915        WGS        298  893305958 PRJNA622674 SAMN14524778 Microbe\\, viral or environmental 340321778
# Center.Name Collection_Date Consent DATASTORE.filetype DATASTORE.provider                     DATASTORE.region
# 1 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   fastq,sra,run.zq         s3,ncbi,gs gs.us-east1,s3.us-east-1,ncbi.public
# 2 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   fastq,run.zq,sra         ncbi,gs,s3 ncbi.public,s3.us-east-1,gs.us-east1
# 3 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   fastq,run.zq,sra         gs,ncbi,s3 ncbi.public,s3.us-east-1,gs.us-east1
# 4 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   sra,run.zq,fastq         s3,gs,ncbi ncbi.public,s3.us-east-1,gs.us-east1
# 5 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   run.zq,fastq,sra         gs,s3,ncbi ncbi.public,gs.us-east1,s3.us-east-1
# 6 NATIONAL INSTITUTE FOR BIOLOGICAL STANDARDS AND CONTROL      2017-01-10  public   run.zq,sra,fastq         ncbi,gs,s3 ncbi.public,gs.us-east1,s3.us-east-1
# Experiment geo_loc_name_country geo_loc_name_country_continent   geo_loc_name  Instrument isolation_source            lat_lon        Library.Name
# 1 SRX8063932       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W 294_294-Sample8_S66
# 2 SRX8063931       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W   306_Sample_16_S16
# 3 SRX8063930       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W   306_Sample_22_S22
# 4 SRX8063929       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W        316_7089_S36
# 5 SRX8063928       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W   306_Sample_24_S24
# 6 SRX8063926       United Kingdom                         Europe United Kingdom NextSeq 500          culture 51.6884 N 0.2409 W   306_Sample_19_S19
# LibraryLayout LibrarySelection LibrarySource                   Organism                              Lifestyle Platform          ReleaseDate
# 1        PAIRED           RANDOM       GENOMIC    Ruminococcus gauvreauii             anaerobic acetate producer ILLUMINA 2020-04-30T00:00:00Z
# 2        PAIRED           RANDOM       GENOMIC     Roseburia intestinalis             thrives on fibre-rich diet ILLUMINA 2020-04-30T00:00:00Z
# 3        PAIRED           RANDOM       GENOMIC          Roseburia hominis             thrives on fibre-rich diet ILLUMINA 2020-04-30T00:00:00Z
# 4        PAIRED           RANDOM       GENOMIC  Prevotella melaninogenica upper respiratory tract, opportunistic ILLUMINA 2020-04-30T00:00:00Z
# 5        PAIRED           RANDOM       GENOMIC            Segatella copri                                   <NA> ILLUMINA 2020-04-30T00:00:00Z
# 6        PAIRED           RANDOM       GENOMIC Parabacteroides distasonis                                   <NA> ILLUMINA 2020-04-30T00:00:00Z
#            create_date version Sample.Name SRA.Study      HOST isolate sample_type  strain sub_species fxn_sum_counts            abbrev     samples variable
# 1 2020-04-07T04:11:00Z       1   DSM_19829 SRP255413 DSM_19829 missing     missing missing        <NA>         981165     R. gauvreauii SRR11487909  Shannon
# 2 2020-04-07T04:32:00Z       1   DSM_14610 SRP255413 DSM_14610 missing     missing missing        <NA>        1027843   R. intestinalis SRR11487910  Shannon
# 3 2020-04-07T04:40:00Z       1   DSM_16839 SRP255413 DSM_16839 missing     missing missing        <NA>        1238611        R. hominis SRR11487911  Shannon
# 4 2020-04-07T04:18:00Z       1    DSM_7089 SRP255413  DSM_7089 missing     missing missing        <NA>        1510358 P. melaninogenica SRR11487912  Shannon
# 5 2020-04-07T04:16:00Z       1   DSM_18205 SRP255413 DSM_18205 missing     missing missing        <NA>        1306958          S. copri SRR11487913  Shannon
# 6 2020-04-07T04:23:00Z       1   DSM_20701 SRP255413 DSM_20701 missing     missing missing        <NA>        1068084     P. distasonis SRR11487915  Shannon
#      value se
# 1 5.177787 NA
# 2 5.339833 NA
# 3 5.290613 NA
# 4 5.215404 NA
# 5 5.308486 NA
# 6 5.401012 NA

# order of samples by Shannon diversity

a_div$data$abbrev[ order(a_div$data$value, decreasing = FALSE)]
# [1] "L. gasseri"             "C. aerofaciens"         "R. gauvreauii"          "B. longum ssp longum"   "P. melaninogenica"      "F. prausnitzii"        
# [7] "B. longum ssp infantis" "B. wexlerae"            "A. hallii"              "A. hadrus"              "R. hominis"             "S. copri"              
# [13] "A. muciniphila"         "A. finegoldii"          "B. thetaiotaomicron"    "R. intestinalis"        "B. uniformis"           "C. butyricum"          
# [19] "P. distasonis"          "E. coli"     

# check: R. gauvreauii SRR11487909 5.177787 ; P. distasonis SRR11487915 5.401012

a_div$data$abbrev <- factor(a_div$data$abbrev, levels = c(a_div$data$abbrev[ order(a_div$data$value, decreasing = FALSE)]), ordered = TRUE )


#set.seed(123)
p <- ggplot(data=a_div$data, aes(x=abbrev, y=value)) +
  #theme_bw()+
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "Shannon diversity") +
  
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(b)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-Alpha-diversity-Shannon-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")



####
#### Assess CPP in selected compounds
####

dat <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-NIBSC_WGS_cultures_PRJNA622674.rds")


## c) Glucose
sel.cpd <- which(df.comp$name == "D-Glucose")
this_var <- "Glucose"

df.comp[sel.cpd, ]
#             id    abbrev      name    form 
# 27    cpd00027     glc-D D-Glucose C6H12O6 
# 24094 cpd26821 D-Glucose D-Glucose C6H12O6   

sel <- which(dat$cpd_id == "cpd00027")

head(dat[sel, ])

temp_dat <- dat[sel, ]
temp_dat$abbrev <- factor(temp_dat$abbrev, levels = c(temp_dat$abbrev[ order(temp_dat$cpd_rel_abun, decreasing = FALSE)]), ordered = TRUE )


p <- ggplot(data=temp_dat, aes(x=abbrev, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "Glucose CPP(%)") +
  
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(c)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-c-Glucose-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")



## d) Cellulose
sel.cpd <- which(df.comp$name == "Cellulose")
this_var <- "Cellulose"

df.comp[sel.cpd, ]
#             id    abbrev      name      form  OC_ratio HC_ratio NC_ratio 
# 11571 cpd11746 Cellulose Cellulose C6H10O5R2 0.8333333 1.666667        0

sel <- which(dat$cpd_id == "cpd11746")

head(dat[sel, ])
temp_dat <- dat[sel, ]
temp_dat$abbrev <- factor(temp_dat$abbrev, levels = c(temp_dat$abbrev[ order(temp_dat$cpd_rel_abun, decreasing = FALSE)]), ordered = TRUE )

p <- ggplot(data=temp_dat, aes(x=abbrev, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "Cellulose CPP(%)") +
  
  theme(
    plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 12.5,unit = "pt"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(d)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-d-Cellulose-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")




## e) CO2 - "Carbon dioxide"                               
sel.cpd <- which(df.comp$name == "CO2")
this_var <- "CO2"

df.comp[sel.cpd, ]
#          id abbrev name
# 11 cpd00011    co2  CO2

sel <- which(dat$cpd_id == "cpd00011")


head(dat[sel, ])
temp_dat <- dat[sel, ]
temp_dat$abbrev <- factor(temp_dat$abbrev, levels = c(temp_dat$abbrev[ order(temp_dat$cpd_rel_abun, decreasing = FALSE)]), ordered = TRUE )

p <- ggplot(data=temp_dat, aes(x=abbrev, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "CO2 CPP(%)") +
  
  theme(
    plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p


grid.text(label = "(e)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-e-CO2-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")



## f) will be heatmap



## g) O2 - "Oxygen"
sel.cpd <- which(df.comp$name == "O2")
this_var <- "O2"

df.comp[sel.cpd, ]
#         id abbrev name form
# 7 cpd00007     o2   O2   O2   

sel <- which(dat$cpd_id == "cpd00007")

head(dat[sel, ])
temp_dat <- dat[sel, ]
temp_dat$abbrev <- factor(temp_dat$abbrev, levels = c(temp_dat$abbrev[ order(temp_dat$cpd_rel_abun, decreasing = FALSE)]), ordered = TRUE )

p <- ggplot(data=temp_dat, aes(x=abbrev, y=cpd_rel_abun)) +
  #theme_bw()+
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "O2 CPP(%)") +
  
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 25,unit = "pt"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p


grid.text(label = "(g)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-g-O2-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")








## h) 
# adenylate energy charge (AEC) indicates the energetic status of soil microorganisms
# the energy status of soilmicroorganisms was evaluated by determining AEC defined as: 
# AEC = (ATP + 0.5 × ADP) / (ATP + ADP + AMP)

sel.cpd <- which(df.comp$name == "ATP")
df.comp[sel.cpd, ]
#.        id abbrev name          form OC_ratio HC_ratio NC_ratio 
# 2 cpd00002    atp  ATP C10H13N5O13P3      1.3      1.3      0.5  
sel <- which(dat$cpd_id == "cpd00002")
head(dat[sel, ])
vals <- list()
vals[["ATP"]] <- dat[sel, ]
vals[["ATP"]]$sample

sel.cpd <- which(df.comp$name == "ADP")
df.comp[sel.cpd, ]
#         id abbrev name          form OC_ratio HC_ratio NC_ratio PC_ratio NP_ratio O_count N_count P_count S_count mass SC_ratio MgC_ratio ZnC_ratio KC_ratio
# 8 cpd00008    adp  ADP C10H13N5O10P2        1      1.3      0.5 
sel <- which(dat$cpd_id == "cpd00008")
head(dat[sel, ])
vals[["ADP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["ADP"]]$sample ) # TRUE

sel.cpd <- which(df.comp$name == "AMP")
df.comp[sel.cpd, ]
#          id abbrev name        form OC_ratio HC_ratio NC_ratio
# 18 cpd00018    amp  AMP C10H12N5O7P      0.7      1.2      0.5 
sel <- which(dat$cpd_id == "cpd00018")
head(dat[sel, ])
vals[["AMP"]] <- dat[sel, ]
identical( vals[["ATP"]]$sample , vals[["AMP"]]$sample ) # TRUE

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
AMP <- vals[["AMP"]]$cpd_rel_abun

AEC <- (ATP + 0.5*ADP) / (ATP + ADP + AMP)

temp <- cbind(dat[sel, ],data.frame(AEC=AEC))
head(temp)

temp$abbrev <- factor(temp$abbrev, levels = c(temp$abbrev[ order(temp$AEC, decreasing = FALSE)]), ordered = TRUE )

p <- ggplot(data=temp, aes(x=abbrev, y=AEC)) +
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  labs(x = NULL, y = "AEC, from CPP(%)") +
  theme(
    plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(h)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-h-AEC-ratio-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")



## i) ATP / ADP

# use vectors for ATP and ADP from above

# calculation
ATP <- vals[["ATP"]]$cpd_rel_abun
ADP <- vals[["ADP"]]$cpd_rel_abun
#AMP <- vals[["AMP"]]$cpd_rel_abun

ATP_ADP_ratio <- ATP/ADP

temp <- cbind(dat[sel, ],data.frame(ATP_ADP_ratio=ATP_ADP_ratio))

temp$abbrev <- factor(temp$abbrev, levels = c(temp$abbrev[ order(temp$ATP_ADP_ratio, decreasing = FALSE)]), ordered = TRUE )


p <- ggplot(data=temp, aes(x=abbrev, y=ATP_ADP_ratio)) +
  theme_classic()+
  geom_point(aes(fill=abbrev), shape = 21, size = 2) + # , alpha = 0.6
  scale_fill_manual(values=cols.bact, name = "Bacteria")+
  
  labs(x = NULL, y = "ATP/ADP, from CPP(%)") +
  theme(
    #plot.margin = margin(t = 5.5,r = 5.5,b = 5.5,l = 18,unit = "pt"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.key.size = unit(0.9,"line"),
    legend.title = element_text(size = rel(0.9)),
    legend.text = element_text(size = rel(0.8)) ,
    axis.text.x  = element_text(angle=30, hjust=1, vjust = 1, size = rel(1))
  )

p

grid.text(label = "(i)", x = unit(0.03, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=14, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-i-ATP_ADP_ratio-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 13, height = 7, units = "cm", res = 600, compression = "lzw",type="cairo")




## f) HEATMAP

dim(dat) # 123840      6
head(dat)
# cpd_id      sample cpd_rel_abun log10_abun                Organism        abbrev
# 1 cpd00035 SRR11487909 2.238630e-01 -0.6500177 Ruminococcus gauvreauii R. gauvreauii
# 2 cpd00117 SRR11487909 2.140593e-01 -0.6694659 Ruminococcus gauvreauii R. gauvreauii
# 3 cpd00041 SRR11487909 2.150706e-01 -0.6674190 Ruminococcus gauvreauii R. gauvreauii
# 4 cpd00320 SRR11487909 1.941327e-05 -4.7119014 Ruminococcus gauvreauii R. gauvreauii
# 5 cpd00504 SRR11487909 5.496867e-02 -1.2598848 Ruminococcus gauvreauii R. gauvreauii
# 6 cpd00516 SRR11487909 9.451592e-02 -1.0244950 Ruminococcus gauvreauii R. gauvreauii

# p<- ggplot(dat, aes(x = sample, y = cpd_id, fill = log10_abun)) + # ggplot(long_df, aes(x = column_name, y = row_id, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_minimal()
# p


#library(pheatmap)

dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ abbrev', value.var = "log10_abun") # use this!
#dat.wide <- reshape2::dcast(dat, formula = 'cpd_id ~ abbrev', value.var = "cpd_rel_abun")

row.names(dat.wide) <- dat.wide[ ,1]
dat.wide <- dat.wide[ ,-1]

# delete compounds with minimal variation
#sel <- which(dat.wide == 0) # 

# calculate row standard deviation?
row_sd <- apply(dat.wide, 1, sd, na.rm = TRUE)
hist(row_sd); summary(row_sd)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.04458 0.74671 1.29639 1.40984 2.12602 3.67952 

dim(dat.wide) # 6192   20
length(row_sd) # 6192
#low_sd_rows <- df[row_sd < 0.5, ]

quantile(row_sd, probs = 0.5)
# 50% 
# 1.296389
sel <- which(row_sd > quantile(row_sd, probs = 0.5) ) # 3095

dat.wide.select <- dat.wide[sel, ]

#pheatmap(as.matrix(dat.wide))
#pheatmap(as.matrix(dat.wide), scale = "row", show_rownames = FALSE)

#pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 6)
pheatmap(as.matrix(dat.wide.select), scale = "row", show_rownames = FALSE, fontsize_col = 10)

grid.text(label = "(f)", x = unit(0.033, "npc") , y = unit(0.97,"npc"), gp=gpar(fontsize=18, fontface="bold") )
#dev.print(tiff, file = paste0(workdir,"/plots/","CPP-f-HEATMAP-REDUCED-SET-facet_grid-v1.tiff"), width = 15, height = 18, units = "cm", res = 500, compression = "lzw",type="cairo")
dev.print(tiff, file = paste0(workdir,"/plots/","CPP-f-HEATMAP-NIBSC_WGS_cultures_PRJNA622674-v8.tiff"), width = 18, height = 22, units = "cm", res = 500, compression = "lzw",type="cairo")


#-------------------------

