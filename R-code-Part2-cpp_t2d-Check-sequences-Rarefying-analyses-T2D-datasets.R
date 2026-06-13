#########################
#
# cpp_t2d: Examining shared trends in compound-associated functional capacities 
# of degraded ecosystem soil microbiomes and type 2 diabetes gut microbiomes
#
# Using compound processing potential (CPP): https://github.com/liddic/cpp_t2d
# Craig Liddicoat | Flinders University, South Australia 
#
# PART 2 - R code to check for sequence counts ~ key signals and re-analysis of rarefied T2D datasets
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


#### Check number of sequences underpinning key signals - 5 soils, 2 T2D datasets - compound groupings - Sugars, Lignin + precursors BCFA-ACPs
#-------------------------


###
### Sunbad Post-mining restoration
###

phy <- readRDS("phy-phyloseq-fxn-sunbad-resto-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-sunbad-resto-v8d.rds")
str(df)
# 'data.frame':	125550 obs. of  6 variables:


## Sunbad Post-mining restoration - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)


df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "Post-mining restoration"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.15257, df = 13, p-value = 0.8811
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4803889  0.5427854
# sample estimates:
#        cor 
# 0.04227791 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Sunbad-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## Sunbad Chronosequence - Sugars
## D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-sunbad-resto-v8d.rds")
str(df) # 'data.frame':	125550 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "Post-mining restoration"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.74477, df = 13, p-value = 0.4697
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3458055  0.6474607
# sample estimates:
#   cor 
# 0.202292 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Sunbad-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Sunbad Chronosequence - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-sunbad-resto-v8d.rds")
str(df) # 'data.frame':	125550 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	45 obs. of  6 variables:
#   $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "mgm4679658.3" "mgm4679658.3" "mgm4679658.3" "mgm4679659.3" ...
# $ cpd_rel_abun: num  0.00199 0.00199 0.00593 0.00209 0.00209 ...
# $ log10_abun  : num  -2.7 -2.7 -2.23 -2.68 -2.68 ...
# $ group       : Ord.factor w/ 5 levels "6"<"12"<"22"<..: 3 3 3 4 4 4 4 4 4 5 ...
# $ ord_group   : int  3 3 3 4 4 4 4 4 4 5 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "Post-mining restoration"

str(res)
# 'data.frame':	15 obs. of  6 variables:
# $ sample      : chr  "mgm4679658.3" "mgm4679659.3" "mgm4679660.3" "mgm4679661.3" ...
# $ sum_rel_abun: num  0.00991 0.01037 0.00912 0.01208 0.00952 ...
# $ group_label : chr  "22" "31" "31" "UM" ...
# $ n_reads     : int  14689920 14381703 14429089 15317615 14581834 17832346 12690709 15318064 12685954 15907230 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "Post-mining restoration" "Post-mining restoration" "Post-mining restoration" "Post-mining restoration" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.4642, df = 13, p-value = 0.6502
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6006634  0.4114881
# sample estimates:
#        cor 
# -0.1276907 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Sunbad-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





###
### AMI
###

phy <- readRDS("phy-phyloseq-fxn-ami-dist-vs-nat-v8c.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-ami-dist-vs-nat-v8c.rds")
str(df) # 'data.frame':	716436 obs. of  5 variables:


## AMI Disturbed vs Natural - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$n_reads_cleaned[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "Disturbed vs Natural"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.67189, df = 82, p-value = 0.5035
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1426636  0.2838864
# sample estimates:
#   cor 
# 0.0739948 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-AMI-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## AMI Disturbed vs Natural - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-ami-dist-vs-nat-v8c.rds")
str(df) # 'data.frame':	716436 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$n_reads_cleaned[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "Disturbed vs Natural"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.12063, df = 82, p-value = 0.9043
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2270667  0.2016511
# sample estimates:
#   cor 
# -0.01331993 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-AMI-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## AMI Disturbed vs Natural - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-ami-dist-vs-nat-v8c.rds")
str(df) # 'data.frame':	716436 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	252 obs. of  5 variables:
# $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "x12465" "x12465" "x12465" "x12469" ...
# $ cpd_rel_abun: num  0.00184 0.00184 0.00549 0.00162 0.00162 ...
# $ log10_abun  : num  -2.74 -2.74 -2.26 -2.79 -2.79 ...
# $ group       : Ord.factor w/ 2 levels "Disturbed"<"Natural": 2 2 2 2 2 2 2 2 2 2 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$n_reads_cleaned[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "Disturbed vs Natural"

str(res)
# 'data.frame':	84 obs. of  6 variables:
# $ sample      : chr  "x12465" "x12469" "x12471" "x12473" ...
# $ sum_rel_abun: num  0.00917 0.00808 0.00863 0.00923 0.00896 ...
# $ group_label : chr  "Natural" "Natural" "Natural" "Natural" ...
# $ n_reads     : int  15589051 15060698 15530417 14967874 16507523 20975240 16895080 16254778 13360505 15197235 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "Disturbed vs Natural" "Disturbed vs Natural" "Disturbed vs Natural" "Disturbed vs Natural" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.37663, df = 82, p-value = 0.7074
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1743933  0.2536906
# sample estimates:
#        cor 
# 0.04155562 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-AMI-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





###
### Mason Prairie restoration
###

phy <- readRDS("phy-phyloseq-fxn-mason-prairie-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-mason-prairie-v8c.rds")
str(df) # 'data.frame':	255810 obs. of  6 variables:


## Mason Prairie restoration - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "Prairie restoration"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 1.0154, df = 28, p-value = 0.3186
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1843288  0.5138378
# sample estimates:
#       cor 
# 0.1884554 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(g)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Mason-BCFA-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## Mason Prairie restoration - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-mason-prairie-v8c.rds")
str(df) # 'data.frame':	255810 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "Prairie restoration"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.75347, df = 28, p-value = 0.4575
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2310315  0.4770141
# sample estimates:
#   cor 
# 0.1409713 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(h)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Mason-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Mason Prairie restoration - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-mason-prairie-v8c.rds")
str(df) # 'data.frame':	255810 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	90 obs. of  6 variables:
# $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "SRR12390916" "SRR12390916" "SRR12390916" "SRR12390917" ...
# $ cpd_rel_abun: num  0.00157 0.00157 0.00467 0.00151 0.00151 ...
# $ log10_abun  : num  -2.81 -2.81 -2.33 -2.82 -2.82 ...
# $ group       : Ord.factor w/ 8 levels "0"<"4"<"5"<"9"<..: 7 7 7 4 4 4 3 3 3 3 ...
# $ ord_group   : num  7 7 7 4 4 4 3 3 3 3 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "Prairie restoration"

str(res)
# 'data.frame':	30 obs. of  6 variables:
# $ sample      : chr  "SRR12390916" "SRR12390917" "SRR12917088" "SRR12917089" ...
# $ sum_rel_abun: num  0.0078 0.00749 0.00621 0.00725 0.007 ...
# $ group_label : chr  "30" "9" "5" "5" ...
# $ n_reads     : int  26451272 33557265 37061571 34628899 44460721 64916268 40708523 66671364 52359741 61396957 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "Prairie restoration" "Prairie restoration" "Prairie restoration" "Prairie restoration" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 1.3674, df = 28, p-value = 0.1824
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1209869  0.5599827
# sample estimates:
#       cor 
# 0.2501874 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(i)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Mason-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



###
### Zuo Plantation succession
###

phy <- readRDS("phy-phyloseq-fxn-zuo_succession-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-zuo-succession-v8c.rds")
str(df) # 'data.frame':	509580 obs. of  6 variables:


## Zuo Plantation succession - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "Plantation succession"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.10391, df = 58, p-value = 0.9176
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2411176  0.2666434
# sample estimates:
#   cor 
# 0.0136424 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(j)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Zuo-BCFA-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## Zuo Plantation succession - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-zuo-succession-v8c.rds")
str(df) # 'data.frame':	509580 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "Plantation succession"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.24066, df = 58, p-value = 0.8107
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2832380  0.2241374
# sample estimates:
#   cor 
# -0.03158488 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(k)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Zuo-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Zuo Plantation succession - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-zuo-succession-v8c.rds")
str(df) # 'data.frame':	509580 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	180 obs. of  6 variables:
#   $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "SRR32132906" "SRR32132906" "SRR32132906" "SRR32132907" ...
# $ cpd_rel_abun: num  0.00113 0.00113 0.00335 0.00114 0.00114 ...
# $ log10_abun  : num  -2.95 -2.95 -2.47 -2.94 -2.94 ...
# $ group       : Ord.factor w/ 4 levels "young"<"half-mature"<..: 2 2 2 2 2 2 2 2 2 2 ...
# $ ord_group   : int  2 2 2 2 2 2 2 2 2 2 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "Plantation succession"

str(res)
# 'data.frame':	60 obs. of  6 variables:
#   $ sample      : chr  "SRR32132906" "SRR32132907" "SRR32132908" "SRR32132909" ...
# $ sum_rel_abun: num  0.00562 0.00564 0.00608 0.00548 0.00437 ...
# $ group_label : chr  "half-mature" "half-mature" "half-mature" "half-mature" ...
# $ n_reads     : int  20155088 19856980 23174145 22823854 18724170 21712139 17944881 20948104 18021696 18266886 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "Plantation succession" "Plantation succession" "Plantation succession" "Plantation succession" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.46833, df = 58, p-value = 0.6413
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1955944  0.3104647
# sample estimates:
#       cor 
# 0.0613788 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(l)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Zuo-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




###
### Wang Veg succession
###


phy <- readRDS("phy-phyloseq-fxn-wang-veg-succession-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-wang-veg-succession-v8c.rds")
str(df) # 'data.frame':	125235 obs. of  6 variables:


## Wang Vegetation succession - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "Vegetation succession"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.052124, df = 13, p-value = 0.9592
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5228463  0.5015215
# sample estimates:
#         cor 
# -0.01445498 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(m)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Wang-BCFA-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## Wang Vegetation succession - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-wang-veg-succession-v8c.rds")
str(df) # 'data.frame':	125235 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "Vegetation succession"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -1.4564, df = 13, p-value = 0.169
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.7440471  0.1704214
# sample estimates:
#   cor 
# -0.3745389 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(n)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Wang-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## Wang Vegetation succession - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-wang-veg-succession-v8c.rds")
str(df) # 'data.frame':	125235 obs. of  6 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	45 obs. of  6 variables:
#   $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "SRR28123688" "SRR28123688" "SRR28123688" "SRR28123689" ...
# $ cpd_rel_abun: num  0.00132 0.00132 0.00392 0.00127 0.00127 ...
# $ log10_abun  : num  -2.88 -2.88 -2.41 -2.9 -2.9 ...
# $ group       : Ord.factor w/ 5 levels "1"<"5"<"20"<"35"<..: 4 4 4 3 3 3 3 3 3 3 ...
# $ ord_group   : int  4 4 4 3 3 3 3 3 3 3 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$clean_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "Vegetation succession"

str(res)
# 'data.frame':	15 obs. of  6 variables:
#   $ sample      : chr  "SRR28123688" "SRR28123689" "SRR28123690" "SRR28123691" ...
# $ sum_rel_abun: num  0.00656 0.00627 0.00643 0.00667 0.00706 ...
# $ group_label : chr  "35" "20" "20" "20" ...
# $ n_reads     : int  12427971 26841936 12206168 11154915 10313095 13883316 12170828 11606332 14413136 15651798 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "Vegetation succession" "Vegetation succession" "Vegetation succession" "Vegetation succession" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.62859, df = 13, p-value = 0.5405
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6286993  0.3733614
# sample estimates:
#        cor 
# -0.1717502 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(o)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-Wang-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




###
### T2D-CHN
###


phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
str(df) # 'data.frame':	591466 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 4.3895, df = 80, p-value = 3.44e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2471905 0.6001903
# sample estimates:
#       cor 
# 0.4405648 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
                      "Pearson cor = ",round(cortest$estimate,3),", P < 0.001"
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(p)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
str(df) # 'data.frame':	591466 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# data:  res$n_reads and res$sum_rel_abun
# t = -3.9537, df = 80, p-value = 0.0001652
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5711890 -0.2052999
# sample estimates:
#   cor 
# -0.4042951 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
                      "Pearson cor = ",round(cortest$estimate,3),", P < 0.001"
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(q)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval.rds")
str(df) # 'data.frame':	591466 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	246 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.45 -8.45 -8.45 -8.45 -8.45 ...
# $ group       : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN"

str(res)
# 'data.frame':	82 obs. of  6 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341586" "SRR341587" ...
# $ sum_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ group_label : chr  "T2D met-" "T2D met-" "T2D met-" "T2D met-" ...
# $ n_reads     : int  12245278 10532029 758416 6751495 5952509 9018931 8257694 12931599 9380015 8014149 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN" "T2D-CHN" "T2D-CHN" "T2D-CHN" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 0.60266, df = 80, p-value = 0.5484
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1519978  0.2801470
# sample estimates:
#        cor 
# 0.06722682 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(r)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-CHN-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")





###
### T2D SWE
###


phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(df) # 'data.frame':	533824 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE"

str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = 1.1547, df = 74, p-value = 0.2519
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.09528012  0.34804964
# sample estimates:
#       cor 
# 0.1330332 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(s)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(df) # 'data.frame':	533824 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE"


str(res)

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -2.1851, df = 74, p-value = 0.03204
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.44685066 -0.02196185
# sample estimates:
#   cor 
# -0.2461984 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(t)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-v6.rds")
str(df) # 533824 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	228 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -8.19 -8.19 -8.19 -8.19 -8.19 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA, n_reads = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  sel.phy <- which(row.names(phy@sam_data) == this_samp)
  res$n_reads[i] <- phy@sam_data$non_host_reads[sel.phy]
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE"

str(res)
# 'data.frame':	76 obs. of  6 variables:
# $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.002731 0.000155 ...
# $ group_label : chr  "T2D met-" "T2D met-" "T2D met-" "Normal" ...
# $ n_reads     : int  5248535 5378909 7906030 6729275 7922852 4516029 6855675 5116818 14235824 19753575 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE" "T2D-SWE" "T2D-SWE" "T2D-SWE" ...

plot(res$n_reads, res$sum_rel_abun)

cortest <- cor.test(x = res$n_reads, y = res$sum_rel_abun)
cortest
# Pearson's product-moment correlation
# data:  res$n_reads and res$sum_rel_abun
# t = -0.74385, df = 74, p-value = 0.4593
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3056678  0.1420657
# sample estimates:
#         cor 
# -0.08614915 

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Pearson cor = ",round(cortest$estimate,3),", P = ",round(cortest$p.value,3)
)

p <- ggplot(data = res, aes(x = n_reads, y = sum_rel_abun) )+
  geom_point()+
  xlab("Number of sequences")+ ylab("Summed CPP (%)")+
  #geom_smooth(method = "lm", alpha = 0.2)+
  geom_smooth(method = "loess", alpha = 0.1)+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "bottom", label = test_result, size = 2.75, lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )
p

grid.text(label = "(u)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Number-sequences-T2D-SWE-Lignin&precursors-v8d.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------



########################## RAREFIED SUBSETS ... 
##########################
########################## T2D-CHN P20
##########################

#### T2D Chinese (CHN) cohort - RERUN subset with even sequences

#### Forslund T2D-CHN - w/ Host-removal - only retain samples with at least 1st quartile (>= 20th percentile) number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-chn.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-chn.rds")

hist(non_host_reads);summary(non_host_reads)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

# only retain samples with at least >= 20th percentile number of sequences

quantile(x = non_host_reads, probs = 0.20)
# 20% 
# 9261724 

length(non_host_reads) # 82

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.20)) # 65

keep_t2d_chn_list_20th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_chn_list])
# SRR413600 SRR341601 SRR341606 SRR413576 SRR341585 SRR341669 SRR413585 SRR413601 SRR413581 SRR413584 SRR413578 SRR341684 SRR413598 SRR341681 SRR413599 SRR341636 SRR413593 
# 9299455   9380015   9622533  10121832  10532029  10557145  10621763  11104392  11172893  11287140  11324956  11378179  11439109  11571973  11587346  11653456  11656012 
# SRR341674 SRR341665 SRR413592 SRR341581 SRR341657 SRR413580 SRR341661 SRR341664 SRR341675 SRR341600 SRR413587 SRR413579 SRR341673 SRR341663 SRR341687 SRR341655 SRR341693 
# 11661126  11771604  12190296  12245278  12461477  12486834  12627861  12773543  12896029  12931599  12985830  13430172  13574397  13796552  13801165  13932846  13980060 
# SRR341676 SRR413610 SRR341713 SRR413575 SRR413626 SRR413618 SRR413617 SRR413625 SRR341670 SRR413608 SRR413621 SRR413614 SRR413670 SRR413606 SRR413615 SRR413603 SRR413637 
# 14455490  15561515  15590289  16621886  16732966  16848086  17288577  17704603  18902323  19002354  19244441  19447073  20330351  20421218  20700045  20778190  21237663 
# SRR413661 SRR413652 SRR413623 SRR413616 SRR413613 SRR413594 SRR413660 SRR413620 SRR413619 SRR413634 SRR413688 SRR413607 SRR413768 SRR413605 
# 21259957  21715990  22081372  22173494  22371270  22470555  23461601  23539875  23569689  23967132  24886198  24958737  26243113  28613098 

writeLines(keep_t2d_chn_list_20th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/keep_t2d_chn_list_20th.txt")

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

sampid <- keep_t2d_chn_list_20th
length(sampid) # 65

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341600"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341601"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341606"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341636"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341585" 
# SRR341600 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341600" 
# SRR341601 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341601" 
# SRR341606 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341606" 
# SRR341636 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_SRR341636" 

length(results_dirs) # 65

sel <- which(names(results_dirs) %in% sampid) # qty 65
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 65

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # FALSE


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 598536      6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                           subsys_L3
# 2 SRR341581 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 3 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine                Glycine Biosynthesis
# 4 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 5 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 6 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 7 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine             Glycine cleavage system
# fxn percent_abun
# 2                                                              Creatinine_amidohydrolase_(EC_3.5.2.10)  0.018162005
# 3                                                           L-threonine_3-dehydrogenase_(EC_1.1.1.103)  0.027243008
# 4                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)  0.015437704
# 5 L-serine_dehydratase,_beta_subunit_(EC_4.3.1.17)_/_L-serine_dehydratase,_alpha_subunit_(EC_4.3.1.17)  0.012713404
# 6                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.012713404
# 7                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.003632401


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 18285    66

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 18285
length(unique(full_fxn_names)) # 18285

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Amino_acid_racemase_RacX" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341600"    "SRR341601"    "SRR341606"    "SRR341636"    "SRR341655"    "SRR341657"    "SRR341661"    "SRR341663"    "SRR341664"   
# [13] "SRR341665"    "SRR341669"    "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"    "SRR341681"    "SRR341684"    "SRR341687"    "SRR341693"    "SRR341713"   
# [25] "SRR413575"    "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413584"    "SRR413585"    "SRR413587"    "SRR413592"    "SRR413593"    "SRR413594"   
# [37] "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"    "SRR413603"    "SRR413605"    "SRR413606"    "SRR413607"    "SRR413608"    "SRR413610"    "SRR413613"    "SRR413614"   
# [49] "SRR413615"    "SRR413616"    "SRR413617"    "SRR413618"    "SRR413619"    "SRR413620"    "SRR413621"    "SRR413623"    "SRR413625"    "SRR413626"    "SRR413634"    "SRR413637"   
# [61] "SRR413652"    "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"    "SRR413768"   

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR413581" "SRR341585" "SRR413584" "SRR413585" "SRR413587"

length(sampid) # 65

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 18285     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 184
length(unique(tax.fxn$subsys_L3)) # 1069
length(unique(tax.fxn$fxn)) # 9721


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 18285     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 18285    65

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18285 taxa and 65 samples ]
# tax_table()   Taxonomy Table:    [ 18285 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341600" "SRR341601" "SRR341606" "SRR341636" "SRR341655" "SRR341657" "SRR341661" "SRR341663" "SRR341664" "SRR341665" "SRR341669" "SRR341670" "SRR341673"
# [16] "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684" "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578" "SRR413579" "SRR413580" "SRR413581" "SRR413584"
# [31] "SRR413585" "SRR413587" "SRR413592" "SRR413593" "SRR413594" "SRR413598" "SRR413599" "SRR413600" "SRR413601" "SRR413603" "SRR413605" "SRR413606" "SRR413607" "SRR413608" "SRR413610"
# [46] "SRR413613" "SRR413614" "SRR413615" "SRR413616" "SRR413617" "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623" "SRR413625" "SRR413626" "SRR413634" "SRR413637" "SRR413652"
# [61] "SRR413660" "SRR413661" "SRR413670" "SRR413688" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse the sample metadata from the non-rarefied phyloseq object

temp <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)

# remove fields that don't pertain to this rarefied data
sel <- which(names(df.samp) %in% c("Bases","total_bases..run.", "non_host_reads", "fxn_sum_counts"))

df.samp <- df.samp[ ,-sel]

# check alignment of names
identical(sample_names(phy), row.names(df.samp)) # TRUE

dim(df.samp) # 65 29


SAMP <- sample_data(df.samp)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18285 taxa and 65 samples ]
# sample_data() Sample Data:       [ 65 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18285 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#       subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Amino_acid_racemase_RacX"      

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   45             20 


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"


saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS")
#phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS")

str(df.samp)
# 'data.frame':	65 obs. of  29 variables:
table( df.samp$gender )
# female   male 
# 31     34 
sel <- which(df.samp$Diagnosis == "T2D metformin-")
table( df.samp$gender[sel] )
# female   male 
# 9     11 
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   59.00   61.00   61.67   63.00   70.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# [1] 9
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   58.00   60.82   69.00   75.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# [1] 11


sel <- which(df.samp$Diagnosis == "ND CTRL")
table( df.samp$gender[sel] )
# female   male 
# 22     23 
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   55.50   56.41   59.75   67.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# [1] 22
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   54.00   57.00   59.04   62.50   74.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# [1] 23

# T2D met- (total n = xx total; females n = x, ages xx-xx; males n = xx, ages xx-xx)
# Normal (total n = xx total; females n = xx, ages xx-xx; males n = xx, ages xx-xx)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 65
NonZeroFxns

mean(NonZeroFxns) # 9208.246
sd(NonZeroFxns) # 2904.958


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 20th percentile)
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
# # For study - Forslund et al T2D-CHN rarefied sequences - 20th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th/working"
# 
# message("\nthis_study <- '-t2d-chn-rarefied-20th-pawsey'")
# this_study <- "-t2d-chn-rarefied-20th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS")
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
#     print(paste0("completed fxn ", f))
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

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 20th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th/working'
# 
# this_study <- '-t2d-chn-rarefied-20th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS')
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
# [1] 18285     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# [[18283]]
# NULL
# 
# [[18284]]
# NULL
# 
# [[18285]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 18285
# [1] "added df 1 of 18285"
# [1] "added df 2 of 18285"
# [1] "added df 3 of 18285"
# ...
# 
# [1] "added df 18283 of 18285"
# [1] "added df 18284 of 18285"
# [1] "added df 18285 of 18285"
# 
# str(df.out)
# 'data.frame':	511744 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 511743      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 10032
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18285 taxa and 65 samples ]
# sample_data() Sample Data:       [ 65 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18285 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.86464
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 10029"
# [1] "completed 10030"
# [1] "completed 10031"
# [1] "completed 10032"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 10032
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341600 SRR341601 SRR341606 SRR341636 SRR341655 SRR341657 
# 100       100       100       100       100       100       100       100 
# SRR341661 SRR341663 SRR341664 SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 
# 100       100       100       100       100       100       100       100 
# SRR341675 SRR341676 SRR341681 SRR341684 SRR341687 SRR341693 SRR341713 SRR413575 
# 100       100       100       100       100       100       100       100 
# SRR413576 SRR413578 SRR413579 SRR413580 SRR413581 SRR413584 SRR413585 SRR413587 
# 100       100       100       100       100       100       100       100 
# SRR413592 SRR413593 SRR413594 SRR413598 SRR413599 SRR413600 SRR413601 SRR413603 
# 100       100       100       100       100       100       100       100 
# SRR413605 SRR413606 SRR413607 SRR413608 SRR413610 SRR413613 SRR413614 SRR413615 
# 100       100       100       100       100       100       100       100 
# SRR413616 SRR413617 SRR413618 SRR413619 SRR413620 SRR413621 SRR413623 SRR413625 
# 100       100       100       100       100       100       100       100 
# SRR413626 SRR413634 SRR413637 SRR413652 SRR413660 SRR413661 SRR413670 SRR413688 
# 100       100       100       100       100       100       100       100 
# SRR413768 
# 100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_20th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 18285    65
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# 
# 
# [[64]]
# NULL
# 
# [[65]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 64"
# [1] "completed 65"
# 
# str(dat)
# 'data.frame':	33263295 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4247.218
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.34181
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7042
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# [[1]]
# NULL
# 
# [[2]]
# NULL
# ...
# 
# 
# 
# [[64]]
# NULL
# 
# [[65]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.12499336
# 4 cpd00117 SRR341581   0.09381302
# 5 cpd00051 SRR341581   0.05902652
# 6 cpd00586 SRR341581   0.00000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 63"
# [1] "completed 64"
# [1] "completed 65"
# 
# str(dat.cpd.collate)
# 'data.frame':	457730 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.125 0.0938 0.059 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4247.218
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.34181
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/cpp_analysis_20th/dat.cpd.collate-all-samps-cpp3d--t2d-chn-rarefied-20th-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	457730 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.125 0.0938 0.059 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000001  0.000118  0.009279  0.001227 11.640452 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -6.164  -3.926    -Inf  -2.911   1.066 


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
# -8.124  -6.164  -3.926  -4.529  -2.911   1.066 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18285 taxa and 65 samples ]
# sample_data() Sample Data:       [ 65 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18285 taxa by 4 taxonomic ranks ]

head(phy@sam_data)
# Sample Data:        [6 samples by 29 sample variables]:
#   Run actual_read_length..run. Age Assay.Type AvgSpotLen  BioProject    BioSample      Bytes center_name..exp. Center.Name Consent DATASTORE.filetype DATASTORE.provider
# SRR341581 SRR341581                      148  59        WGS        148 PRJNA422434 SAMN00715131 1170792822               BGI         BGI  public          fastq,sra                 s3
# SRR341585 SRR341585                      148  60        WGS        148 PRJNA422434 SAMN00715135 1200370883               BGI         BGI  public          sra,fastq                 s3
# SRR341600 SRR341600                      148  70        WGS        148 PRJNA422434 SAMN00715150 1254879513               BGI         BGI  public          sra,fastq                 s3
# SRR341601 SRR341601                      148  68        WGS        148 PRJNA422434 SAMN00715151  948467380               BGI         BGI  public          fastq,sra                 s3
# SRR341606 SRR341606                      148  63        WGS        148 PRJNA422434 SAMN00715156 1024883659               BGI         BGI  public          fastq,sra                 s3
# SRR341636 SRR341636                      148  51        WGS        148 PRJNA422434 SAMN00715186 1075830990               BGI         BGI  public          fastq,sra                 s3
# DATASTORE.region Experiment gender                  Instrument     Library.Name LibraryLayout LibrarySelection LibrarySource NATION             Organism Platform          ReleaseDate
# SRR341581     s3.us-east-1  SRX095662 female Illumina Genome Analyzer II HGMlijMCFDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# SRR341585     s3.us-east-1  SRX095666 female Illumina Genome Analyzer II HGMlijMDGDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# SRR341600     s3.us-east-1  SRX095681 female Illumina Genome Analyzer II HGMlijMCZDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# SRR341601     s3.us-east-1  SRX095682 female Illumina Genome Analyzer II HGMlijMBYDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# SRR341606     s3.us-east-1  SRX095687 female Illumina Genome Analyzer II HGMlijMDODFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# SRR341636     s3.us-east-1  SRX095717 female Illumina Genome Analyzer II HGMlijMBRDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z
# run..run. Sample.Name SRA.Study      Diagnosis
# SRR341581 FC615J5AAXX  bgi-DLF001 SRP008047 T2D metformin-
#   SRR341585 FC61B1KAAXX  bgi-DLF005 SRP008047 T2D metformin-
#   SRR341600 FC615J5AAXX  bgi-DOF007 SRP008047 T2D metformin-
#   SRR341601 FC61B1DAAXX  bgi-DOF009 SRP008047 T2D metformin-
#   SRR341606 FC61B1KAAXX  bgi-DOF014 SRP008047 T2D metformin-
#   SRR341636 FC61B1DAAXX  bgi-NOF002 SRP008047        ND CTRL

samp <- as(phy@sam_data,"data.frame")
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   
samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group_label[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"  
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-20th-qty65-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-20th-qty65-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	457730 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.125 0.0938 0.059 ...
# $ log10_abun  : num  -8.124 -8.124 -0.903 -1.028 -1.229 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 7042
7042*65 # 457730


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 SRR341581   0.00000000 -8.1236827    T2D met-
#   2 cpd02175 SRR341581   0.00000000 -8.1236827    T2D met-
#   3 cpd00035 SRR341581   0.12499336 -0.9031131    T2D met-
#   4 cpd00117 SRR341581   0.09381302 -1.0277369    T2D met-
#   5 cpd00051 SRR341581   0.05902652 -1.2289528    T2D met-
#   6 cpd00586 SRR341581   0.00000000 -8.1236827    T2D met-

dim(data_in) # 457730      5

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

mean(no_compounds) # 5707.662
sd(no_compounds) #  1059.63

mean(sample_sum_relabun) # 65.34181
sd(sample_sum_relabun) # 3.910869

length(unique(data_in$cpd_id)) # 7042

#-------------------------



# all p < 0.05
#### Forslund T2D-CHN - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 20th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty65-EVEN-seqs-20th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-20th-qty65-v8e.rds")
str(df) # 'data.frame':	457730 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2340 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.12 -8.12 -8.12 -8.12 -8.12 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	65 obs. of  5 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341600" "SRR341601" ...
# $ sum_rel_abun: num  0.00 0.00 0.00 0.00 1.67e-05 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 2 1 1 1 1 ...
# $ cpd_group   : chr  "BCFA-ACPs" "BCFA-ACPs" "BCFA-ACPs" "BCFA-ACPs" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" ...

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 20
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 43

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 213, p-value = 0.0003871
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,5))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  expand_limits(y = 1.2*max(res$sum_rel_abun))+
  #ylim( min(res$sum_rel_abun), 0.015 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-20th-qty65-v8e.rds")
str(df) # 'data.frame':	457730 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	325 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0.4798 0.0886 0.0871 0.2175 0.0902 ...
# $ log10_abun  : num  -0.319 -1.053 -1.06 -0.663 -1.045 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	65 obs. of  5 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341600" "SRR341601" ...
# $ sum_rel_abun: num  0.963 1.109 0.756 1.378 1.41 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 2 1 1 1 1 ...
# $ cpd_group   : chr  "Sugars" "Sugars" "Sugars" "Sugars" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" ...

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 642, p-value = 0.002897
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 1.5 )+
  expand_limits(y = 1.15*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-20th-qty65-v8e.rds")
str(df) # 457730 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	195 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.12 -8.12 -8.12 -8.12 -8.12 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	65 obs. of  5 variables:
# $ sample      : chr  "SRR341581" "SRR341585" "SRR341600" "SRR341601" ...
# $ sum_rel_abun: num  0 0 0 0 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 2 1 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" "T2D-CHN Rarefied (P20)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -4.888  -4.357    -Inf  -3.791  -2.815 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.824  -4.888  -4.357  -4.399  -3.791  -2.815 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 309, p-value = 0.02277
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.5 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-CHN-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------



##########################
########################## T2D-CHN P15
##########################
##########################

#### T2D Chinese (CHN) cohort - RERUN subset with even sequences

#### Forslund T2D-CHN - w/ Host-removal - only retain samples with at least >= 15th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-chn.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-chn.rds")

hist(non_host_reads);summary(non_host_reads)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

# only retain samples with at least 1st quartile (>= 15th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.15)
# 15% 
# 8050681 

length(non_host_reads) # 82

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.15)) # 69

keep_t2d_chn_list_15th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_chn_list_15th])
# SRR341599 SRR341589 SRR413758 SRR341652 SRR413600 SRR341601 SRR341606 SRR413576 SRR341585 SRR341669 SRR413585 SRR413601 SRR413581 SRR413584 SRR413578 SRR341684 SRR413598 
# 8257694   9018931   9094862   9252291   9299455   9380015   9622533  10121832  10532029  10557145  10621763  11104392  11172893  11287140  11324956  11378179  11439109 
# SRR341681 SRR413599 SRR341636 SRR413593 SRR341674 SRR341665 SRR413592 SRR341581 SRR341657 SRR413580 SRR341661 SRR341664 SRR341675 SRR341600 SRR413587 SRR413579 SRR341673 
# 11571973  11587346  11653456  11656012  11661126  11771604  12190296  12245278  12461477  12486834  12627861  12773543  12896029  12931599  12985830  13430172  13574397 
# SRR341663 SRR341687 SRR341655 SRR341693 SRR341676 SRR413610 SRR341713 SRR413575 SRR413626 SRR413618 SRR413617 SRR413625 SRR341670 SRR413608 SRR413621 SRR413614 SRR413670 
# 13796552  13801165  13932846  13980060  14455490  15561515  15590289  16621886  16732966  16848086  17288577  17704603  18902323  19002354  19244441  19447073  20330351 
# SRR413606 SRR413615 SRR413603 SRR413637 SRR413661 SRR413652 SRR413623 SRR413616 SRR413613 SRR413594 SRR413660 SRR413620 SRR413619 SRR413634 SRR413688 SRR413607 SRR413768 
# 20421218  20700045  20778190  21237663  21259957  21715990  22081372  22173494  22371270  22470555  23461601  23539875  23569689  23967132  24886198  24958737  26243113 
# SRR413605 
# 28613098 

writeLines(keep_t2d_chn_list_15th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/keep_t2d_chn_list_15th.txt")

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

sampid <- keep_t2d_chn_list_15th
length(sampid) # 69

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341589"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341599"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341600"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341601"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341585" 
# SRR341589 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341589" 
# SRR341599 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341599" 
# SRR341600 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341600" 
# SRR341601 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_SRR341601" 

length(results_dirs) # 69

sel <- which(names(results_dirs) %in% sampid) # qty 69
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 69

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 599665      6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                           subsys_L3
# 2 SRR341581 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 3 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine                Glycine Biosynthesis
# 4 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 5 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 6 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 7 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine             Glycine cleavage system
# fxn percent_abun
# 2                                                              Creatinine_amidohydrolase_(EC_3.5.2.10)  0.020872469
# 3                                                           L-threonine_3-dehydrogenase_(EC_1.1.1.103)  0.031308704
# 4                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)  0.017741599
# 5 L-serine_dehydratase,_beta_subunit_(EC_4.3.1.17)_/_L-serine_dehydratase,_alpha_subunit_(EC_4.3.1.17)  0.014610728
# 6                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.014610728
# 7                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.004174494


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 18160    70

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 18160
length(unique(full_fxn_names)) # 18160

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Amino_acid_racemase_RacX" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341589"    "SRR341599"    "SRR341600"    "SRR341601"    "SRR341606"    "SRR341636"    "SRR341652"    "SRR341655"    "SRR341657"   
# [13] "SRR341661"    "SRR341663"    "SRR341664"    "SRR341665"    "SRR341669"    "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"    "SRR341681"    "SRR341684"   
# [25] "SRR341687"    "SRR341693"    "SRR341713"    "SRR413575"    "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413584"    "SRR413585"    "SRR413587"   
# [37] "SRR413592"    "SRR413593"    "SRR413594"    "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"    "SRR413603"    "SRR413605"    "SRR413606"    "SRR413607"    "SRR413608"   
# [49] "SRR413610"    "SRR413613"    "SRR413614"    "SRR413615"    "SRR413616"    "SRR413617"    "SRR413618"    "SRR413619"    "SRR413620"    "SRR413621"    "SRR413623"    "SRR413625"   
# [61] "SRR413626"    "SRR413634"    "SRR413637"    "SRR413652"    "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"    "SRR413758"    "SRR413768" 

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR413581" "SRR341585" "SRR413584" "SRR413585" "SRR341589"

length(sampid) # 69

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 18160     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 182
length(unique(tax.fxn$subsys_L3)) # 1063
length(unique(tax.fxn$fxn)) # 9640


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 18160     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 18160    69

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18160 taxa and 69 samples ]
# tax_table()   Taxonomy Table:    [ 18160 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341589" "SRR341599" "SRR341600" "SRR341601" "SRR341606" "SRR341636" "SRR341652" "SRR341655" "SRR341657" "SRR341661" "SRR341663" "SRR341664" "SRR341665" "SRR341669"
# [17] "SRR341670" "SRR341673" "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684" "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578" "SRR413579" "SRR413580" "SRR413581"
# [33] "SRR413584" "SRR413585" "SRR413587" "SRR413592" "SRR413593" "SRR413594" "SRR413598" "SRR413599" "SRR413600" "SRR413601" "SRR413603" "SRR413605" "SRR413606" "SRR413607" "SRR413608" "SRR413610"
# [49] "SRR413613" "SRR413614" "SRR413615" "SRR413616" "SRR413617" "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623" "SRR413625" "SRR413626" "SRR413634" "SRR413637" "SRR413652" "SRR413660"
# [65] "SRR413661" "SRR413670" "SRR413688" "SRR413758" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse the sample metadata from the non-rarefied phyloseq object

temp <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)

# remove fields that don't pertain to this rarefied data
sel <- which(names(df.samp) %in% c("Bases","total_bases..run.", "non_host_reads", "fxn_sum_counts"))

df.samp <- df.samp[ ,-sel]

# check alignment of names
identical(sample_names(phy), row.names(df.samp)) # TRUE

dim(df.samp) # 69 29


SAMP <- sample_data(df.samp)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18160 taxa and 69 samples ]
# sample_data() Sample Data:       [ 69 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18160 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Amino_acid_racemase_RacX"     

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   47             22 


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"


saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS")
#phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS")

str(df.samp)
# 'data.frame':	69 obs. of  29 variables:
table( df.samp$gender )
# female   male 
# 34     35 
sel <- which(df.samp$Diagnosis == "T2D metformin-")
table( df.samp$gender[sel] )
# female   male 
#   11     11 
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   57.50   60.00   60.18   63.00   70.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# [1] 11
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   58.00   60.82   69.00   75.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# [1] 11


sel <- which(df.samp$Diagnosis == "ND CTRL")
table( df.samp$gender[sel] )
# female   male 
# 23     24 
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   56.00   56.39   59.50   67.00
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# [1] 23
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   53.75   56.50   58.75   62.25   74.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# [1] 24

# T2D met- (total n = .. total; females n = ., ages ..-..; males n = .., ages ..-..)
# Normal (total n = .. total; females n = .., ages ..-..; males n = .., ages ..-..)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 69
NonZeroFxns

mean(NonZeroFxns) # 8690.797
sd(NonZeroFxns) # 3256.726


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 15th percentile)
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
# # For study - Forslund et al T2D-CHN rarefied sequences - 15th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th/working"
# 
# message("\nthis_study <- '-t2d-chn-rarefied-15th-pawsey'")
# this_study <- "-t2d-chn-rarefied-15th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS")
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
#     print(paste0("completed fxn ", f))
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

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 15th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th/working'
# 
# this_study <- '-t2d-chn-rarefied-15th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS')
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
# [1] 18160     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# 
# [[18158]]
# NULL
# 
# [[18159]]
# NULL
# 
# [[18160]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 18160
# [1] "added df 1 of 18160"
# [1] "added df 2 of 18160"
# [1] "added df 3 of 18160"
# ...
# 
# [1] "added df 18157 of 18160"
# [1] "added df 18158 of 18160"
# [1] "added df 18159 of 18160"
# [1] "added df 18160 of 18160"
# 
# str(df.out)
# 'data.frame':	505344 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 505343      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9968
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18160 taxa and 69 samples ]
# sample_data() Sample Data:       [ 69 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18160 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.88987
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 9965"
# [1] "completed 9966"
# [1] "completed 9967"
# [1] "completed 9968"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9968
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341589 SRR341599 SRR341600 SRR341601 SRR341606 SRR341636 
# 100       100       100       100       100       100       100       100 
# SRR341652 SRR341655 SRR341657 SRR341661 SRR341663 SRR341664 SRR341665 SRR341669 
# 100       100       100       100       100       100       100       100 
# SRR341670 SRR341673 SRR341674 SRR341675 SRR341676 SRR341681 SRR341684 SRR341687 
# 100       100       100       100       100       100       100       100 
# SRR341693 SRR341713 SRR413575 SRR413576 SRR413578 SRR413579 SRR413580 SRR413581 
# 100       100       100       100       100       100       100       100 
# SRR413584 SRR413585 SRR413587 SRR413592 SRR413593 SRR413594 SRR413598 SRR413599 
# 100       100       100       100       100       100       100       100 
# SRR413600 SRR413601 SRR413603 SRR413605 SRR413606 SRR413607 SRR413608 SRR413610 
# 100       100       100       100       100       100       100       100 
# SRR413613 SRR413614 SRR413615 SRR413616 SRR413617 SRR413618 SRR413619 SRR413620 
# 100       100       100       100       100       100       100       100 
# SRR413621 SRR413623 SRR413625 SRR413626 SRR413634 SRR413637 SRR413652 SRR413660 
# 100       100       100       100       100       100       100       100 
# SRR413661 SRR413670 SRR413688 SRR413758 SRR413768 
# 100       100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_15th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 18160    69
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[68]]
# NULL
# 
# [[69]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 67"
# [1] "completed 68"
# [1] "completed 69"
# 
# str(dat)
# 'data.frame':	34868667 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4530.734
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.66281
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7034
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
# ...
# 
# 
# 
# [[68]]
# NULL
# 
# [[69]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.12170224
# 4 cpd00117 SRR341581   0.08879049
# 5 cpd00051 SRR341581   0.05218117
# 6 cpd00586 SRR341581   0.00000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 67"
# [1] "completed 68"
# [1] "completed 69"
# 
# str(dat.cpd.collate)
# 'data.frame':	485346 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1217 0.0888 0.0522 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4530.734
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.66281
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/cpp_analysis_15th/dat.cpd.collate-all-samps-cpp3d--t2d-chn-rarefied-15th-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	485346 obs. of  3 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1217 0.0888 0.0522 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000109  0.009335  0.001182 11.663503 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -6.513  -3.961    -Inf  -2.927   1.067 


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
# -8.083  -6.513  -3.961  -4.606  -2.927   1.067 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18160 taxa and 69 samples ]
# sample_data() Sample Data:       [ 69 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18160 taxa by 4 taxonomic ranks ]

head(phy@sam_data)
# Sample Data:        [6 samples by 29 sample variables]:
#   Run actual_read_length..run. Age Assay.Type AvgSpotLen  BioProject    BioSample      Bytes center_name..exp. Center.Name Consent DATASTORE.filetype DATASTORE.provider DATASTORE.region
# SRR341581 SRR341581                      148  59        WGS        148 PRJNA422434 SAMN00715131 1170792822               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# SRR341585 SRR341585                      148  60        WGS        148 PRJNA422434 SAMN00715135 1200370883               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341589 SRR341589                      148  51        WGS        148 PRJNA422434 SAMN00715139 1010953818               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341599 SRR341599                      148  56        WGS        148 PRJNA422434 SAMN00715149  945051953               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# SRR341600 SRR341600                      148  70        WGS        148 PRJNA422434 SAMN00715150 1254879513               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341601 SRR341601                      148  68        WGS        148 PRJNA422434 SAMN00715151  948467380               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# Experiment gender                  Instrument     Library.Name LibraryLayout LibrarySelection LibrarySource NATION             Organism Platform          ReleaseDate   run..run. Sample.Name
# SRR341581  SRX095662 female Illumina Genome Analyzer II HGMlijMCFDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DLF001
# SRR341585  SRX095666 female Illumina Genome Analyzer II HGMlijMDGDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF005
# SRR341589  SRX095670 female Illumina Genome Analyzer II HGMlijMDSDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF010
# SRR341599  SRX095680 female Illumina Genome Analyzer II HGMlijMDIDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DOF002
# SRR341600  SRX095681 female Illumina Genome Analyzer II HGMlijMCZDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DOF007
# SRR341601  SRX095682 female Illumina Genome Analyzer II HGMlijMBYDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1DAAXX  bgi-DOF009
# SRA.Study      Diagnosis
# SRR341581 SRP008047 T2D metformin-
#   SRR341585 SRP008047 T2D metformin-
#   SRR341589 SRP008047 T2D metformin-
#   SRR341599 SRP008047 T2D metformin-
#   SRR341600 SRP008047 T2D metformin-
#   SRR341601 SRP008047 T2D metformin-

samp <- as(phy@sam_data,"data.frame")
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   
samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group_label[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"  
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-15th-qty69-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-15th-qty69-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	485346 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1217 0.0888 0.0522 ...
# $ log10_abun  : num  -8.083 -8.083 -0.915 -1.052 -1.282 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 7060
7034*69 # 485346


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 SRR341581   0.00000000 -8.0833604    T2D met-
# 2 cpd02175 SRR341581   0.00000000 -8.0833604    T2D met-
# 3 cpd00035 SRR341581   0.12170224 -0.9147014    T2D met-
# 4 cpd00117 SRR341581   0.08879049 -1.0516335    T2D met-
# 5 cpd00051 SRR341581   0.05218117 -1.2824862    T2D met-
# 6 cpd00586 SRR341581   0.00000000 -8.0833604    T2D met-

dim(data_in) # 485346      5

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

mean(no_compounds) # 5528.246
sd(no_compounds) #  1271.285

mean(sample_sum_relabun) # 65.66281
sd(sample_sum_relabun) # 4.200844

length(unique(data_in$cpd_id)) # 7034

#-------------------------


# all p < 0.05
#### Forslund T2D-CHN - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 15th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty69-EVEN-seqs-15th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-15th-qty69-v8e.rds")
str(df) # 'data.frame':	485346 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2484 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.08 -8.08 -8.08 -8.08 -8.08 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] #
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 245, p-value = 0.000236
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,5))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.015 )+
  expand_limits(y = 1.2*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-15th-qty69-v8e.rds")
str(df) # 'data.frame':	485346 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	345 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0.4426 0.0948 0.0931 0.2307 0.0967 ...
# $ log10_abun  : num  -0.354 -1.023 -1.031 -0.637 -1.015 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 747, p-value = 0.001333
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.58 )+
  expand_limits(y = 1.12*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-15th-qty69-v8e.rds")
str(df) # 485346 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	207 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.08 -8.08 -8.08 -8.08 -8.08 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	69 obs. of  5 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341589" "SRR341599" ...
# $ sum_rel_abun: num  0 0 0 0 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 2 2 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P15)" "T2D-CHN Rarefied (P15)" "T2D-CHN Rarefied (P15)" "T2D-CHN Rarefied (P15)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -5.091  -4.371    -Inf  -3.840  -2.811 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.764  -5.091  -4.371  -4.440  -3.840  -2.811 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 339.5, p-value = 0.01116
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.5 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-CHN-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------




##########################
########################## T2D-CHN P10
##########################
##########################

#### T2D Chinese (CHN) cohort - RERUN subset with even sequences

#### Forslund T2D-CHN - w/ Host-removal - only retain samples with at least >= 10th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-chn.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-chn.rds")

hist(non_host_reads);summary(non_host_reads)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

# only retain samples with at least 1st quartile (>= 10th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.10)
# 10% 
# 7356437 

length(non_host_reads) # 82

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.10)) # 73

keep_t2d_chn_list_10th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_chn_list_10th])
# SRR413642 SRR413597 SRR341660 SRR341602 SRR341599 SRR341589 SRR413758 SRR341652 SRR413600 SRR341601 SRR341606 SRR413576 SRR341585 SRR341669 SRR413585 SRR413601 SRR413581 SRR413584 
# 7577408   7752136   7763661   8014149   8257694   9018931   9094862   9252291   9299455   9380015   9622533  10121832  10532029  10557145  10621763  11104392  11172893  11287140 
# SRR413578 SRR341684 SRR413598 SRR341681 SRR413599 SRR341636 SRR413593 SRR341674 SRR341665 SRR413592 SRR341581 SRR341657 SRR413580 SRR341661 SRR341664 SRR341675 SRR341600 SRR413587 
# 11324956  11378179  11439109  11571973  11587346  11653456  11656012  11661126  11771604  12190296  12245278  12461477  12486834  12627861  12773543  12896029  12931599  12985830 
# SRR413579 SRR341673 SRR341663 SRR341687 SRR341655 SRR341693 SRR341676 SRR413610 SRR341713 SRR413575 SRR413626 SRR413618 SRR413617 SRR413625 SRR341670 SRR413608 SRR413621 SRR413614 
# 13430172  13574397  13796552  13801165  13932846  13980060  14455490  15561515  15590289  16621886  16732966  16848086  17288577  17704603  18902323  19002354  19244441  19447073 
# SRR413670 SRR413606 SRR413615 SRR413603 SRR413637 SRR413661 SRR413652 SRR413623 SRR413616 SRR413613 SRR413594 SRR413660 SRR413620 SRR413619 SRR413634 SRR413688 SRR413607 SRR413768 
# 20330351  20421218  20700045  20778190  21237663  21259957  21715990  22081372  22173494  22371270  22470555  23461601  23539875  23569689  23967132  24886198  24958737  26243113 
# SRR413605 
# 28613098 

writeLines(keep_t2d_chn_list_10th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/keep_t2d_chn_list_10th.txt")

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

sampid <- keep_t2d_chn_list_10th
length(sampid) # 73

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341589"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341599"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341600"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341601"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341585" 
# SRR341589 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341589" 
# SRR341599 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341599" 
# SRR341600 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341600" 
# SRR341601 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_SRR341601" 

length(results_dirs) # 73

sel <- which(names(results_dirs) %in% sampid) # qty 73
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 73

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 620337      6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                           subsys_L3
# 2 SRR341581 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation
# 3 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine                Glycine Biosynthesis
# 4 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 5 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 6 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine      Glycine and Serine Utilization
# 7 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine             Glycine cleavage system
# fxn percent_abun
# 2                                                              Creatinine_amidohydrolase_(EC_3.5.2.10)  0.022810219
# 3                                                           L-threonine_3-dehydrogenase_(EC_1.1.1.103)  0.034215328
# 4                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)  0.009124088
# 5 L-serine_dehydratase,_beta_subunit_(EC_4.3.1.17)_/_L-serine_dehydratase,_alpha_subunit_(EC_4.3.1.17)  0.015967153
# 6                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.015967153
# 7                                                                   L-serine_dehydratase_(EC_4.3.1.17)  0.004562044


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 18089    74

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 18089
length(unique(full_fxn_names)) # 18089

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Amino_acid_racemase_RacX" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341589"    "SRR341599"    "SRR341600"    "SRR341601"    "SRR341602"    "SRR341606"    "SRR341636"    "SRR341652"    "SRR341655"   
# [13] "SRR341657"    "SRR341660"    "SRR341661"    "SRR341663"    "SRR341664"    "SRR341665"    "SRR341669"    "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"   
# [25] "SRR341681"    "SRR341684"    "SRR341687"    "SRR341693"    "SRR341713"    "SRR413575"    "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413584"   
# [37] "SRR413585"    "SRR413587"    "SRR413592"    "SRR413593"    "SRR413594"    "SRR413597"    "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"    "SRR413603"    "SRR413605"   
# [49] "SRR413606"    "SRR413607"    "SRR413608"    "SRR413610"    "SRR413613"    "SRR413614"    "SRR413615"    "SRR413616"    "SRR413617"    "SRR413618"    "SRR413619"    "SRR413620"   
# [61] "SRR413621"    "SRR413623"    "SRR413625"    "SRR413626"    "SRR413634"    "SRR413637"    "SRR413642"    "SRR413652"    "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"   
# [73] "SRR413758"    "SRR413768"

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR413581" "SRR341585" "SRR413584" "SRR413585" "SRR341589"

length(sampid) # 73

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 18089     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 183
length(unique(tax.fxn$subsys_L3)) # 1063
length(unique(tax.fxn$fxn)) # 9610


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 18089     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 18089    73

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18089 taxa and 73 samples ]
# tax_table()   Taxonomy Table:    [ 18089 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341589" "SRR341599" "SRR341600" "SRR341601" "SRR341602" "SRR341606" "SRR341636" "SRR341652" "SRR341655" "SRR341657" "SRR341660" "SRR341661" "SRR341663" "SRR341664"
# [17] "SRR341665" "SRR341669" "SRR341670" "SRR341673" "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684" "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578" "SRR413579"
# [33] "SRR413580" "SRR413581" "SRR413584" "SRR413585" "SRR413587" "SRR413592" "SRR413593" "SRR413594" "SRR413597" "SRR413598" "SRR413599" "SRR413600" "SRR413601" "SRR413603" "SRR413605" "SRR413606"
# [49] "SRR413607" "SRR413608" "SRR413610" "SRR413613" "SRR413614" "SRR413615" "SRR413616" "SRR413617" "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623" "SRR413625" "SRR413626" "SRR413634"
# [65] "SRR413637" "SRR413642" "SRR413652" "SRR413660" "SRR413661" "SRR413670" "SRR413688" "SRR413758" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse the sample metadata from the non-rarefied phyloseq object

temp <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)

# remove fields that don't pertain to this rarefied data
sel <- which(names(df.samp) %in% c("Bases","total_bases..run.", "non_host_reads", "fxn_sum_counts"))

df.samp <- df.samp[ ,-sel]

# check alignment of names
identical(sample_names(phy), row.names(df.samp)) # TRUE

dim(df.samp) # 73 29


SAMP <- sample_data(df.samp)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18089 taxa and 73 samples ]
# sample_data() Sample Data:       [ 73 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18089 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Amino_acid_racemase_RacX"           

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   49             24 


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"


saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS")
#phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS")

str(df.samp)
# 'data.frame':	73 obs. of  29 variables:
table( df.samp$gender )
# female   male 
# 35     38 
sel <- which(df.samp$Diagnosis == "T2D metformin-")
table( df.samp$gender[sel] )
# female   male 
# 12     12 
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   58.25   60.50   60.33   63.00   70.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# [1] 12
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   57.00   60.33   68.50   75.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# [1] 12


sel <- which(df.samp$Diagnosis == "ND CTRL")
table( df.samp$gender[sel] )
# female   male 
# 23     26 
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   56.00   56.39   59.50   67.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# [1] 23
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   53.25   56.00   58.38   61.50   74.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# [1] 26

# T2D met- (total n = .. total; females n = ., ages ..-..; males n = .., ages ..-..)
# Normal (total n = .. total; females n = .., ages ..-..; males n = .., ages ..-..)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 73
NonZeroFxns

mean(NonZeroFxns) # 8497.767
sd(NonZeroFxns) # 3277.732


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 10th percentile)
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
# # For study - Forslund et al T2D-CHN rarefied sequences - 10th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th/working"
# 
# message("\nthis_study <- '-t2d-chn-rarefied-10th-pawsey'")
# this_study <- "-t2d-chn-rarefied-10th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS")
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
#     print(paste0("completed fxn ", f))
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

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 10th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th/working'
# 
# this_study <- '-t2d-chn-rarefied-10th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS')
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
# [1] 18089     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# [[18088]]
# NULL
# 
# [[18089]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 18089
# [1] "added df 1 of 18089"
# [1] "added df 2 of 18089"
# [1] "added df 3 of 18089"
# ...
# 
# [1] "added df 18087 of 18089"
# [1] "added df 18088 of 18089"
# [1] "added df 18089 of 18089"
# 
# str(df.out)
# 'data.frame':	505549 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 505548      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9925
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18089 taxa and 73 samples ]
# sample_data() Sample Data:       [ 73 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18089 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.8676
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 9923"
# [1] "completed 9924"
# [1] "completed 9925"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9925
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341589 SRR341599 SRR341600 SRR341601 SRR341602 SRR341606 
# 100       100       100       100       100       100       100       100 
# SRR341636 SRR341652 SRR341655 SRR341657 SRR341660 SRR341661 SRR341663 SRR341664 
# 100       100       100       100       100       100       100       100 
# SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 SRR341675 SRR341676 SRR341681 
# 100       100       100       100       100       100       100       100 
# SRR341684 SRR341687 SRR341693 SRR341713 SRR413575 SRR413576 SRR413578 SRR413579 
# 100       100       100       100       100       100       100       100 
# SRR413580 SRR413581 SRR413584 SRR413585 SRR413587 SRR413592 SRR413593 SRR413594 
# 100       100       100       100       100       100       100       100 
# SRR413597 SRR413598 SRR413599 SRR413600 SRR413601 SRR413603 SRR413605 SRR413606 
# 100       100       100       100       100       100       100       100 
# SRR413607 SRR413608 SRR413610 SRR413613 SRR413614 SRR413615 SRR413616 SRR413617 
# 100       100       100       100       100       100       100       100 
# SRR413618 SRR413619 SRR413620 SRR413621 SRR413623 SRR413625 SRR413626 SRR413634 
# 100       100       100       100       100       100       100       100 
# SRR413637 SRR413642 SRR413652 SRR413660 SRR413661 SRR413670 SRR413688 SRR413758 
# 100       100       100       100       100       100       100       100 
# SRR413768 
# 100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_10th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 18089    73
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[72]]
# NULL
# 
# [[73]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 71"
# [1] "completed 72"
# [1] "completed 73"
# 
# str(dat)
# 'data.frame':	36905004 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4801.174
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.76951
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7019
# 
# ### 3) collate_compounds within each sample
# 
# # # # # # # # # # #
# [[1]]
# NULL
# 
# [[2]]
# NULL
# ...
# 
# 
# 
# [[72]]
# NULL
# 
# [[73]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.11575981
# 4 cpd00117 SRR341581   0.09010390
# 5 cpd00051 SRR341581   0.05132299
# 6 cpd00586 SRR341581   0.00000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 71"
# [1] "completed 72"
# [1] "completed 73"
# 
# str(dat.cpd.collate)
# 'data.frame':	512387 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1158 0.0901 0.0513 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4801.174
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.76951
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/cpp_analysis_10th/dat.cpd.collate-all-samps-cpp3d--t2d-chn-rarefied-10th-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	512387 obs. of  3 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1158 0.0901 0.0513 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000109  0.009370  0.001190 11.705165 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -6.557  -3.962    -Inf  -2.924   1.068 


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
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -8.044  -6.557  -3.962  -4.605  -2.924   1.068 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 18089 taxa and 73 samples ]
# sample_data() Sample Data:       [ 73 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 18089 taxa by 4 taxonomic ranks ]

head(phy@sam_data)
# includes field for Diagnosis
# Sample Data:        [6 samples by 29 sample variables]:
#   Run actual_read_length..run. Age Assay.Type AvgSpotLen  BioProject    BioSample      Bytes center_name..exp. Center.Name Consent DATASTORE.filetype DATASTORE.provider DATASTORE.region
# SRR341581 SRR341581                      148  59        WGS        148 PRJNA422434 SAMN00715131 1170792822               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# SRR341585 SRR341585                      148  60        WGS        148 PRJNA422434 SAMN00715135 1200370883               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341589 SRR341589                      148  51        WGS        148 PRJNA422434 SAMN00715139 1010953818               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341599 SRR341599                      148  56        WGS        148 PRJNA422434 SAMN00715149  945051953               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# SRR341600 SRR341600                      148  70        WGS        148 PRJNA422434 SAMN00715150 1254879513               BGI         BGI  public          sra,fastq                 s3     s3.us-east-1
# SRR341601 SRR341601                      148  68        WGS        148 PRJNA422434 SAMN00715151  948467380               BGI         BGI  public          fastq,sra                 s3     s3.us-east-1
# Experiment gender                  Instrument     Library.Name LibraryLayout LibrarySelection LibrarySource NATION             Organism Platform          ReleaseDate   run..run. Sample.Name
# SRR341581  SRX095662 female Illumina Genome Analyzer II HGMlijMCFDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DLF001
# SRR341585  SRX095666 female Illumina Genome Analyzer II HGMlijMDGDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF005
# SRR341589  SRX095670 female Illumina Genome Analyzer II HGMlijMDSDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF010
# SRR341599  SRX095680 female Illumina Genome Analyzer II HGMlijMDIDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DOF002
# SRR341600  SRX095681 female Illumina Genome Analyzer II HGMlijMCZDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DOF007
# SRR341601  SRX095682 female Illumina Genome Analyzer II HGMlijMBYDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA 2012-09-05T00:00:00Z FC61B1DAAXX  bgi-DOF009
# SRA.Study      Diagnosis
# SRR341581 SRP008047 T2D metformin-
#   SRR341585 SRP008047 T2D metformin-
#   SRR341589 SRP008047 T2D metformin-
#   SRR341599 SRP008047 T2D metformin-
#   SRR341600 SRP008047 T2D metformin-
#   SRR341601 SRP008047 T2D metformin-

samp <- as(phy@sam_data,"data.frame")
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   
samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group_label[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"  
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-10th-qty73-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-10th-qty73-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	512387 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1158 0.0901 0.0513 ...
# $ log10_abun  : num  -8.044 -8.044 -0.936 -1.045 -1.29 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 7019
7019*73 # 512387


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 SRR341581   0.00000000 -8.0437793    T2D met-
# 2 cpd02175 SRR341581   0.00000000 -8.0437793    T2D met-
# 3 cpd00035 SRR341581   0.11575981 -0.9364422    T2D met-
# 4 cpd00117 SRR341581   0.09010390 -1.0452564    T2D met-
# 5 cpd00051 SRR341581   0.05132299 -1.2896880    T2D met-
# 6 cpd00586 SRR341581   0.00000000 -8.0437793    T2D met-

dim(data_in) # 512387      5

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

mean(no_compounds) # 5493.589
sd(no_compounds) #  1275.574

mean(sample_sum_relabun) # 65.76951
sd(sample_sum_relabun) # 4.170144

length(unique(data_in$cpd_id)) # 7019

#-------------------------


# all p < 0.05
#### Forslund T2D-CHN - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 10th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty73-EVEN-seqs-10th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-10th-qty73-v8e.rds")
str(df) # 'data.frame':	512387 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2628 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.04 -8.04 -8.04 -8.04 -8.04 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 24
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 49

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 254, p-value = 4.492e-05
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,7))


p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  ylim( min(res$sum_rel_abun), 0.015 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-10th-qty73-v8e.rds")
str(df) # 'data.frame':	512387 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	365 obs. of  5 variables:
# $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0.4504 0.0989 0.097 0.2392 0.1009 ...
# $ log10_abun  : num  -0.346 -1.005 -1.013 -0.621 -0.996 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 861, p-value = 0.0005527
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,5))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.58 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-10th-qty73-v8e.rds")
str(df) # 512387 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	219 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -8.04 -8.04 -8.04 -8.04 -8.04 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	73 obs. of  5 variables:
# $ sample      : chr  "SRR341581" "SRR341585" "SRR341589" "SRR341599" ...
# $ sum_rel_abun: num  0 0 0 0 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 2 2 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P10)" "T2D-CHN Rarefied (P10)" "T2D-CHN Rarefied (P10)" "T2D-CHN Rarefied (P10)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -5.227  -4.412    -Inf  -3.871  -2.814 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.456  -5.227  -4.412  -4.618  -3.871  -2.814 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 417, p-value = 0.02226
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.5 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-CHN-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------



##########################
########################## T2D-CHN P5
##########################
##########################

#### T2D Chinese (CHN) cohort - RERUN subset with even sequences

#### Forslund T2D-CHN - w/ Host-removal - only retain samples with at least >= 5th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-chn.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-chn.rds")

hist(non_host_reads);summary(non_host_reads)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

sum(non_host_reads) # 1159922740 = 1,159,922,740
length(non_host_reads) # 82

# only retain samples with at least >= 5th percentile number of sequences

quantile(x = non_host_reads, probs = 0.05)
# 5% 
# 5962972 

length(non_host_reads) # 82

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.05)) # 77

keep_t2d_chn_list_5th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_chn_list_5th])
# SRR413582 SRR341587 SRR341654 SRR341604 SRR413642 SRR413597 SRR341660 SRR341602 SRR341599 SRR341589 SRR413758 SRR341652 SRR413600 SRR341601 SRR341606 SRR413576 SRR341585 
# 6161771   6751495   7204681   7331885   7577408   7752136   7763661   8014149   8257694   9018931   9094862   9252291   9299455   9380015   9622533  10121832  10532029 
# SRR341669 SRR413585 SRR413601 SRR413581 SRR413584 SRR413578 SRR341684 SRR413598 SRR341681 SRR413599 SRR341636 SRR413593 SRR341674 SRR341665 SRR413592 SRR341581 SRR341657 
# 10557145  10621763  11104392  11172893  11287140  11324956  11378179  11439109  11571973  11587346  11653456  11656012  11661126  11771604  12190296  12245278  12461477 
# SRR413580 SRR341661 SRR341664 SRR341675 SRR341600 SRR413587 SRR413579 SRR341673 SRR341663 SRR341687 SRR341655 SRR341693 SRR341676 SRR413610 SRR341713 SRR413575 SRR413626 
# 12486834  12627861  12773543  12896029  12931599  12985830  13430172  13574397  13796552  13801165  13932846  13980060  14455490  15561515  15590289  16621886  16732966 
# SRR413618 SRR413617 SRR413625 SRR341670 SRR413608 SRR413621 SRR413614 SRR413670 SRR413606 SRR413615 SRR413603 SRR413637 SRR413661 SRR413652 SRR413623 SRR413616 SRR413613 
# 16848086  17288577  17704603  18902323  19002354  19244441  19447073  20330351  20421218  20700045  20778190  21237663  21259957  21715990  22081372  22173494  22371270 
# SRR413594 SRR413660 SRR413620 SRR413619 SRR413634 SRR413688 SRR413607 SRR413768 SRR413605 
# 22470555  23461601  23539875  23569689  23967132  24886198  24958737  26243113  28613098 

writeLines(keep_t2d_chn_list_5th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/keep_t2d_chn_list_5th.txt")

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

sampid <- keep_t2d_chn_list_5th
length(sampid) # 77

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341587"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341589"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341599"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341600"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341585" 
# SRR341587 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341587" 
# SRR341589 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341589" 
# SRR341599 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341599" 
# SRR341600 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_SRR341600" 

length(results_dirs) # 77

sel <- which(names(results_dirs) %in% sampid) # qty 77
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 77

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 624664      6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                      subsys_L3
# 2 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine           Glycine Biosynthesis
# 3 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine Glycine and Serine Utilization
# 4 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine Glycine and Serine Utilization
# 5 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine Glycine and Serine Utilization
# 6 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine        Glycine cleavage system
# 7 SRR341581 Amino Acids and Derivatives Alanine, serine, and glycine            Serine Biosynthesis
# fxn percent_abun
# 2                                                           L-threonine_3-dehydrogenase_(EC_1.1.1.103)   0.02829655
# 3                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)   0.01131862
# 4 L-serine_dehydratase,_beta_subunit_(EC_4.3.1.17)_/_L-serine_dehydratase,_alpha_subunit_(EC_4.3.1.17)   0.01980758
# 5                                                                   L-serine_dehydratase_(EC_4.3.1.17)   0.01980758
# 6                                                                   L-serine_dehydratase_(EC_4.3.1.17)   0.00565931
# 7                                                     D-3-phosphoglycerate_dehydrogenase_(EC_1.1.1.95)   0.01131862


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 17903    78

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 17903
length(unique(full_fxn_names)) # 17903

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylaconitate_cis-trans_isomerase" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Amino_acid_racemase_RacX" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341587"    "SRR341589"    "SRR341599"    "SRR341600"    "SRR341601"    "SRR341602"    "SRR341604"    "SRR341606"   
# [12] "SRR341636"    "SRR341652"    "SRR341654"    "SRR341655"    "SRR341657"    "SRR341660"    "SRR341661"    "SRR341663"    "SRR341664"    "SRR341665"    "SRR341669"   
# [23] "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"    "SRR341681"    "SRR341684"    "SRR341687"    "SRR341693"    "SRR341713"    "SRR413575"   
# [34] "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413582"    "SRR413584"    "SRR413585"    "SRR413587"    "SRR413592"    "SRR413593"   
# [45] "SRR413594"    "SRR413597"    "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"    "SRR413603"    "SRR413605"    "SRR413606"    "SRR413607"    "SRR413608"   
# [56] "SRR413610"    "SRR413613"    "SRR413614"    "SRR413615"    "SRR413616"    "SRR413617"    "SRR413618"    "SRR413619"    "SRR413620"    "SRR413621"    "SRR413623"   
# [67] "SRR413625"    "SRR413626"    "SRR413634"    "SRR413637"    "SRR413642"    "SRR413652"    "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"    "SRR413758"   
# [78] "SRR413768" 

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR413581" "SRR341585" "SRR413582" "SRR341587" "SRR413584"

length(sampid) # 77

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 17903     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 183
length(unique(tax.fxn$subsys_L3)) # 1056
length(unique(tax.fxn$fxn)) # 9523


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 17903     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 17903    77

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17903 taxa and 77 samples ]
# tax_table()   Taxonomy Table:    [ 17903 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341587" "SRR341589" "SRR341599" "SRR341600" "SRR341601" "SRR341602" "SRR341604" "SRR341606" "SRR341636" "SRR341652" "SRR341654" "SRR341655"
# [15] "SRR341657" "SRR341660" "SRR341661" "SRR341663" "SRR341664" "SRR341665" "SRR341669" "SRR341670" "SRR341673" "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684"
# [29] "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578" "SRR413579" "SRR413580" "SRR413581" "SRR413582" "SRR413584" "SRR413585" "SRR413587" "SRR413592"
# [43] "SRR413593" "SRR413594" "SRR413597" "SRR413598" "SRR413599" "SRR413600" "SRR413601" "SRR413603" "SRR413605" "SRR413606" "SRR413607" "SRR413608" "SRR413610" "SRR413613"
# [57] "SRR413614" "SRR413615" "SRR413616" "SRR413617" "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623" "SRR413625" "SRR413626" "SRR413634" "SRR413637" "SRR413642"
# [71] "SRR413652" "SRR413660" "SRR413661" "SRR413670" "SRR413688" "SRR413758" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse the sample metadata from the non-rarefied phyloseq object

temp <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)

# remove fields that don't pertain to this rarefied data
sel <- which(names(df.samp) %in% c("Bases","total_bases..run.", "non_host_reads", "fxn_sum_counts"))

df.samp <- df.samp[ ,-sel]

# check alignment of names
identical(sample_names(phy), row.names(df.samp)) # TRUE

dim(df.samp) # 77 29


SAMP <- sample_data(df.samp)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17903 taxa and 77 samples ]
# sample_data() Sample Data:       [ 77 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 17903 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylaconitate_cis-trans_isomerase"        
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Amino_acid_racemase_RacX"     

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#   50             27 


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"


saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS")
#phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS")

str(df.samp)
# 'data.frame':	77 obs. of  29 variables:
table( df.samp$gender )
# female   male 
# 39     38 
sel <- which(df.samp$Diagnosis == "T2D metformin-")
table( df.samp$gender[sel] )
# female   male 
# 15     12 
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.0    57.5    61.0    60.0    62.5    70.0 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# [1] 15
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   57.00   60.33   68.50   75.00  
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# [1] 12


sel <- which(df.samp$Diagnosis == "ND CTRL")
table( df.samp$gender[sel] )
# female   male 
# 24     26 
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   56.00   56.58   60.00   67.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# [1] 24
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   53.25   56.00   58.38   61.50   74.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# [1] 26

# T2D met- (total n = .. total; females n = ., ages ..-..; males n = .., ages ..-..)
# Normal (total n = .. total; females n = .., ages ..-..; males n = .., ages ..-..)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 77
NonZeroFxns

mean(NonZeroFxns) # 8112.519
sd(NonZeroFxns) # 3388.785


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 5th percentile)
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
# # For study - Forslund et al T2D-CHN rarefied sequences - 5th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th/working"
# 
# message("\nthis_study <- '-t2d-chn-rarefied-5th-pawsey'")
# this_study <- "-t2d-chn-rarefied-5th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS")
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
#     print(paste0("completed fxn ", f))
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

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 5th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th/working'
# 
# this_study <- '-t2d-chn-rarefied-5th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS')
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
# [1] 17903     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[17901]]
# NULL
# 
# [[17902]]
# NULL
# 
# [[17903]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 17903
# [1] "added df 1 of 17903"
# [1] "added df 2 of 17903"
# [1] "added df 3 of 17903"
# ...
# 
# [1] "added df 17900 of 17903"
# [1] "added df 17901 of 17903"
# [1] "added df 17902 of 17903"
# [1] "added df 17903 of 17903"
# 
# str(df.out)
# 'data.frame':	500570 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_2" "fxn_2" "fxn_3" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 500569      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9827
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17903 taxa and 77 samples ]
# sample_data() Sample Data:       [ 77 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 17903 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 54.89024
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 9825"
# [1] "completed 9826"
# [1] "completed 9827"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9827
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341587 SRR341589 SRR341599 SRR341600 SRR341601 SRR341602 
# 100       100       100       100       100       100       100       100 
# SRR341604 SRR341606 SRR341636 SRR341652 SRR341654 SRR341655 SRR341657 SRR341660 
# 100       100       100       100       100       100       100       100 
# SRR341661 SRR341663 SRR341664 SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 
# 100       100       100       100       100       100       100       100 
# SRR341675 SRR341676 SRR341681 SRR341684 SRR341687 SRR341693 SRR341713 SRR413575 
# 100       100       100       100       100       100       100       100 
# SRR413576 SRR413578 SRR413579 SRR413580 SRR413581 SRR413582 SRR413584 SRR413585 
# 100       100       100       100       100       100       100       100 
# SRR413587 SRR413592 SRR413593 SRR413594 SRR413597 SRR413598 SRR413599 SRR413600 
# 100       100       100       100       100       100       100       100 
# SRR413601 SRR413603 SRR413605 SRR413606 SRR413607 SRR413608 SRR413610 SRR413613 
# 100       100       100       100       100       100       100       100 
# SRR413614 SRR413615 SRR413616 SRR413617 SRR413618 SRR413619 SRR413620 SRR413621 
# 100       100       100       100       100       100       100       100 
# SRR413623 SRR413625 SRR413626 SRR413634 SRR413637 SRR413642 SRR413652 SRR413660 
# 100       100       100       100       100       100       100       100 
# SRR413661 SRR413670 SRR413688 SRR413758 SRR413768 
# 100       100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_5th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 17903    77
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# 
# [[76]]
# NULL
# 
# [[77]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_2 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_3 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 75"
# [1] "completed 76"
# [1] "completed 77"
# 
# str(dat)
# 'data.frame':	38543813 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_2" "fxn_2" "fxn_3" "fxn_3" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5059.74
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.7109
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 7076
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
# ...
# 
# 
# 
# [[76]]
# NULL
# 
# [[77]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.11533640
# 4 cpd00117 SRR341581   0.06877136
# 5 cpd00051 SRR341581   0.07074137
# 6 cpd00586 SRR341581   0.00000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 75"
# [1] "completed 76"
# [1] "completed 77"
# 
# str(dat.cpd.collate)
# 'data.frame':	544852 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1153 0.0688 0.0707 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5059.74
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.7109
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/cpp_analysis_5th/dat.cpd.collate-all-samps-cpp3d--t2d-chn-rarefied-5th-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	544852 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1153 0.0688 0.0707 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000101  0.009286  0.001136 11.709863 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -7.174  -3.997    -Inf  -2.945   1.069 


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
# -7.953  -7.174  -3.997  -4.669  -2.945   1.069 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 17903 taxa and 77 samples ]
# sample_data() Sample Data:       [ 77 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 17903 taxa by 4 taxonomic ranks ]

head(phy@sam_data)
# includes field for Diagnosis
# Sample Data:        [6 samples by 29 sample variables]:
#   Run actual_read_length..run. Age Assay.Type AvgSpotLen  BioProject    BioSample      Bytes center_name..exp. Center.Name Consent DATASTORE.filetype DATASTORE.provider
# SRR341581 SRR341581                      148  59        WGS        148 PRJNA422434 SAMN00715131 1170792822               BGI         BGI  public          fastq,sra                 s3
# SRR341585 SRR341585                      148  60        WGS        148 PRJNA422434 SAMN00715135 1200370883               BGI         BGI  public          sra,fastq                 s3
# SRR341587 SRR341587                      148  62        WGS        148 PRJNA422434 SAMN00715137 1479828698               BGI         BGI  public          fastq,sra                 s3
# SRR341589 SRR341589                      148  51        WGS        148 PRJNA422434 SAMN00715139 1010953818               BGI         BGI  public          sra,fastq                 s3
# SRR341599 SRR341599                      148  56        WGS        148 PRJNA422434 SAMN00715149  945051953               BGI         BGI  public          fastq,sra                 s3
# SRR341600 SRR341600                      148  70        WGS        148 PRJNA422434 SAMN00715150 1254879513               BGI         BGI  public          sra,fastq                 s3
# DATASTORE.region Experiment gender                  Instrument     Library.Name LibraryLayout LibrarySelection LibrarySource NATION             Organism Platform
# SRR341581     s3.us-east-1  SRX095662 female Illumina Genome Analyzer II HGMlijMCFDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341585     s3.us-east-1  SRX095666 female Illumina Genome Analyzer II HGMlijMDGDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341587     s3.us-east-1  SRX095668 female Illumina Genome Analyzer II HGMlijMAPDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341589     s3.us-east-1  SRX095670 female Illumina Genome Analyzer II HGMlijMDSDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341599     s3.us-east-1  SRX095680 female Illumina Genome Analyzer II HGMlijMDIDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341600     s3.us-east-1  SRX095681 female Illumina Genome Analyzer II HGMlijMCZDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# ReleaseDate   run..run. Sample.Name SRA.Study      Diagnosis
# SRR341581 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DLF001 SRP008047 T2D metformin-
#   SRR341585 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF005 SRP008047 T2D metformin-
#   SRR341587 2012-09-05T00:00:00Z FC61B81AAXX  bgi-DLF007 SRP008047 T2D metformin-
#   SRR341589 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF010 SRP008047 T2D metformin-
#   SRR341599 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DOF002 SRP008047 T2D metformin-
#   SRR341600 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DOF007 SRP008047 T2D metformin-

samp <- as(phy@sam_data,"data.frame")
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   
samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group_label[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"  
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-5th-qty77-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-5th-qty77-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	544852 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1153 0.0688 0.0707 ...
# $ log10_abun  : num  -7.953 -7.953 -0.938 -1.163 -1.15 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 7076
7076*77 # 544852


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 SRR341581   0.00000000 -8.0437793    T2D met-
# 2 cpd02175 SRR341581   0.00000000 -8.0437793    T2D met-
# 3 cpd00035 SRR341581   0.11575981 -0.9364422    T2D met-
# 4 cpd00117 SRR341581   0.09010390 -1.0452564    T2D met-
# 5 cpd00051 SRR341581   0.05132299 -1.2896880    T2D met-
# 6 cpd00586 SRR341581   0.00000000 -8.0437793    T2D met-

dim(data_in) # 544852      5

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

mean(no_compounds) # 5340.571
sd(no_compounds) #  1403.449

mean(sample_sum_relabun) # 65.7109
sd(sample_sum_relabun) # 4.251117

length(unique(data_in$cpd_id)) # 7076

#-------------------------

# 2 of 3 have P < 0.05
#### Forslund T2D-CHN - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 5th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty77-EVEN-seqs-5th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-5th-qty77-v8e.rds")
str(df) # 'data.frame':	544852 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2772 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -7.95 -7.95 -7.95 -7.95 -7.95 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 27
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 50

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 280, p-value = 1.264e-05
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,7))


p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.015 )+
  expand_limits(y = 1.2*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-5th-qty77-v8e.rds")
str(df) # 'data.frame':	544852 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	385 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0.4177 0.0952 0.0928 0.2315 0.1004 ...
# $ log10_abun  : num  -0.379 -1.021 -1.033 -0.635 -0.998 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 986, p-value = 0.0004588
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,5))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 1.6 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-5th-qty77-v8e.rds")
str(df) # 544852 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	231 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -7.95 -7.95 -7.95 -7.95 -7.95 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	77 obs. of  5 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341587" "SRR341589" ...
# $ sum_rel_abun: num  0 0 0 0 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN Rarefied (P5)" "T2D-CHN Rarefied (P5)" "T2D-CHN Rarefied (P5)" "T2D-CHN Rarefied (P5)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf  -5.162  -4.386    -Inf  -3.803  -2.703 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.633  -5.162  -4.386  -4.401  -3.803  -2.703 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 527.5, p-value = 0.05748
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-CHN-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------



##########################
########################## T2D-CHN MIN
##########################
##########################

#### T2D Chinese (CHN) cohort - RERUN subset with even sequences

#### Forslund T2D-CHN - w/ Host-removal - keep all samples but rarefy to minimum library size
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-chn.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-chn.rds")

hist(non_host_reads);summary(non_host_reads)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 758416 10224381 12700702 14145399 19183919 28613098 

sum(non_host_reads) # 1159922740 = 1,159,922,740
length(non_host_reads) # 82


min(non_host_reads) # 758416


keep_t2d_chn_list_min <- names(non_host_reads)


writeLines(keep_t2d_chn_list_min, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/keep_t2d_chn_list_min.txt")

#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (minimum library size)
#-------------------------

sampid <- keep_t2d_chn_list_min
length(sampid) # 82

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min"

list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# # don't keep 1st two 
# ( results_dirs <- list.dirs(superfocus_out_dir)[-c(1,2)] )

# # don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341581"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341585"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341586"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341587"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341588"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341589"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# SRR341581 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341581" 
# SRR341585 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341585" 
# SRR341586 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341586" 
# SRR341587 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341587" 
# SRR341588 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341588" 
# SRR341589 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_SRR341589" 

length(results_dirs) # 82

sel <- which(names(results_dirs) %in% sampid) # qty 82
#results_dirs <- results_dirs[sel]

length( which(names(results_dirs) %in% sampid)) # 82

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "R1.good.fastq..$", x = names(tab))
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 432269      6
head(sfx.long)
# sampleID                   subsys_L1                                                         subsys_L2                                                   subsys_L3
# 2 SRR341581 Amino Acids and Derivatives                                  Arginine; urea cycle, polyamines                          Arginine and Ornithine Degradation
# 3 SRR341581 Amino Acids and Derivatives                              Aromatic amino acids and derivatives        Aromatic amino acid interconversions with aryl acids
# 4 SRR341581 Amino Acids and Derivatives                              Aromatic amino acids and derivatives        Aromatic amino acid interconversions with aryl acids
# 5 SRR341581 Amino Acids and Derivatives                              Aromatic amino acids and derivatives        Aromatic amino acid interconversions with aryl acids
# 6 SRR341581 Amino Acids and Derivatives                              Aromatic amino acids and derivatives                                        Tryptophan synthesis
# 7 SRR341581 Amino Acids and Derivatives Glutamine, glutamate, aspartate, asparagine; ammonia assimilation Glutamine, Glutamate, Aspartate and Asparagine Biosynthesis
# fxn percent_abun
# 2                                                   Arginase_(EC_3.5.3.1)    0.2141328
# 3                                                       Coenzyme_A_ligase    0.1070664
# 4                           Phenylacetate-coenzyme_A_ligase_(EC_6.2.1.30)    0.1070664
# 5 Phenylacetate-coenzyme_A_ligase_(EC_6.2.1.30);_Coenzyme_F390_synthetase    0.2141328
# 6                            Tryptophan_synthase_beta_chain_(EC_4.2.1.20)    0.2141328
# 7              Asparagine_synthetase_[glutamine-hydrolyzing]_(EC_6.3.5.4)    0.2141328


sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)


## translate from long to wide format

names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 15270    83

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 15270
length(unique(full_fxn_names)) # 15270

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Aspartate_racemase_(EC_5.1.1.13)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Diaminopimelate_epimerase_(EC_5.1.1.7)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)


head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "SRR341581"    "SRR341585"    "SRR341586"    "SRR341587"    "SRR341588"    "SRR341589"    "SRR341599"    "SRR341600"    "SRR341601"    "SRR341602"   
# [12] "SRR341604"    "SRR341606"    "SRR341609"    "SRR341636"    "SRR341645"    "SRR341646"    "SRR341652"    "SRR341654"    "SRR341655"    "SRR341657"    "SRR341660"   
# [23] "SRR341661"    "SRR341663"    "SRR341664"    "SRR341665"    "SRR341669"    "SRR341670"    "SRR341673"    "SRR341674"    "SRR341675"    "SRR341676"    "SRR341681"   
# [34] "SRR341684"    "SRR341687"    "SRR341693"    "SRR341713"    "SRR413575"    "SRR413576"    "SRR413578"    "SRR413579"    "SRR413580"    "SRR413581"    "SRR413582"   
# [45] "SRR413584"    "SRR413585"    "SRR413587"    "SRR413592"    "SRR413593"    "SRR413594"    "SRR413597"    "SRR413598"    "SRR413599"    "SRR413600"    "SRR413601"   
# [56] "SRR413603"    "SRR413605"    "SRR413606"    "SRR413607"    "SRR413608"    "SRR413610"    "SRR413613"    "SRR413614"    "SRR413615"    "SRR413616"    "SRR413617"   
# [67] "SRR413618"    "SRR413619"    "SRR413620"    "SRR413621"    "SRR413623"    "SRR413625"    "SRR413626"    "SRR413634"    "SRR413637"    "SRR413642"    "SRR413652"   
# [78] "SRR413660"    "SRR413661"    "SRR413670"    "SRR413688"    "SRR413758"    "SRR413768" 

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)


head(sampid)
# "SRR341581" "SRR413581" "SRR341585" "SRR413582" "SRR341586" "SRR341587"

length(sampid) # 82

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

# identical(names(sfx.wide) , as.character(gsub(pattern = "-",replacement = "_",x = sampid))) # FALSE
# length( names(sfx.wide) %in% as.character(gsub(pattern = "-",replacement = "_",x = sampid)) ) # 113 - i.e. matching but order different

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 15270     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 176
length(unique(tax.fxn$subsys_L3)) # 985
length(unique(tax.fxn$fxn)) # 8263


#-------------------------

#### Forslund-CHN-T2D - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (minimum library size)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 15270     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 15270    82

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15270 taxa and 82 samples ]
# tax_table()   Taxonomy Table:    [ 15270 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "SRR341581" "SRR341585" "SRR341586" "SRR341587" "SRR341588" "SRR341589" "SRR341599" "SRR341600" "SRR341601" "SRR341602" "SRR341604" "SRR341606" "SRR341609" "SRR341636"
# [15] "SRR341645" "SRR341646" "SRR341652" "SRR341654" "SRR341655" "SRR341657" "SRR341660" "SRR341661" "SRR341663" "SRR341664" "SRR341665" "SRR341669" "SRR341670" "SRR341673"
# [29] "SRR341674" "SRR341675" "SRR341676" "SRR341681" "SRR341684" "SRR341687" "SRR341693" "SRR341713" "SRR413575" "SRR413576" "SRR413578" "SRR413579" "SRR413580" "SRR413581"
# [43] "SRR413582" "SRR413584" "SRR413585" "SRR413587" "SRR413592" "SRR413593" "SRR413594" "SRR413597" "SRR413598" "SRR413599" "SRR413600" "SRR413601" "SRR413603" "SRR413605"
# [57] "SRR413606" "SRR413607" "SRR413608" "SRR413610" "SRR413613" "SRR413614" "SRR413615" "SRR413616" "SRR413617" "SRR413618" "SRR413619" "SRR413620" "SRR413621" "SRR413623"
# [71] "SRR413625" "SRR413626" "SRR413634" "SRR413637" "SRR413642" "SRR413652" "SRR413660" "SRR413661" "SRR413670" "SRR413688" "SRR413758" "SRR413768"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse the sample metadata from the non-rarefied phyloseq object

temp <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)

# remove fields that don't pertain to this rarefied data
sel <- which(names(df.samp) %in% c("Bases","total_bases..run.", "non_host_reads", "fxn_sum_counts"))

df.samp <- df.samp[ ,-sel]

# check alignment of names
identical(sample_names(phy), row.names(df.samp)) # TRUE

dim(df.samp) # 82 29


SAMP <- sample_data(df.samp)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15270 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 15270 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                            
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"      
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"   
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Aspartate_racemase_(EC_5.1.1.13)"             
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Diaminopimelate_epimerase_(EC_5.1.1.7)"   

table(phy@sam_data$Diagnosis)
# ND CTRL T2D metformin- 
#      52             30


getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/cpp3d/modelling/R"


saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS")
#phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS")

str(df.samp)
# 'data.frame':	77 obs. of  29 variables:
table( df.samp$gender )
# female   male 
# 41     41 
sel <- which(df.samp$Diagnosis == "T2D metformin-")
table( df.samp$gender[sel] )
# female   male 
# 17     13 
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   56.00   61.00   60.29   63.00   71.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "female")] )
# [1] 17
summary( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   58.00   61.15   70.00   75.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "T2D metformin-" & df.samp$gender == "male")] )
# [1] 13


sel <- which(df.samp$Diagnosis == "ND CTRL")
table( df.samp$gender[sel] )
# female   male 
# 24     28 
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 51.00   53.00   56.00   56.58   60.00   67.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "female")] )
# [1] 24
summary( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 52.00   53.75   56.00   58.71   62.25   74.00 
length( df.samp$Age[ which(df.samp$Diagnosis == "ND CTRL" & df.samp$gender == "male")] )
# [1] 28

# T2D met- (total n = .. total; females n = ., ages ..-..; males n = .., ages ..-..)
# Normal (total n = .. total; females n = .., ages ..-..; males n = .., ages ..-..)


# get stats??
head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 82
NonZeroFxns

mean(NonZeroFxns) # 5271.573
sd(NonZeroFxns) # 2863.603


#-------------------------

#### Forslund-CHN-T2D - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (minimum library size)
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
# # For study - Forslund et al T2D-CHN rarefied sequences - minimum library size
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min/working"
# 
# message("\nthis_study <- '-t2d-chn-rarefied-min-pawsey'")
# this_study <- "-t2d-chn-rarefied-min-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS")
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
#     print(paste0("completed fxn ", f))
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

#### Forslund-CHN-T2D - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (minimum library size)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min/working'
# 
# this_study <- '-t2d-chn-rarefied-min-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS')
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
# [1] 15270     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[15268]]
# NULL
# 
# [[15269]]
# NULL
# 
# [[15270]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 15270
# [1] "added df 1 of 15270"
# [1] "added df 2 of 15270"
# [1] "added df 3 of 15270"
# ...
# 
# 
# [1] "added df 15268 of 15270"
# [1] "added df 15269 of 15270"
# [1] "added df 15270 of 15270"
# 
# str(df.out)
# 'data.frame':	409980 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_2" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_1  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_1  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 409979      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 8451
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15270 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 15270 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 55.34381
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 8449"
# [1] "completed 8450"
# [1] "completed 8451"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 8451
# 
# sample_sums(phy)
# SRR341581 SRR341585 SRR341586 SRR341587 SRR341588 SRR341589 SRR341599 SRR341600 
# 100       100       100       100       100       100       100       100 
# SRR341601 SRR341602 SRR341604 SRR341606 SRR341609 SRR341636 SRR341645 SRR341646 
# 100       100       100       100       100       100       100       100 
# SRR341652 SRR341654 SRR341655 SRR341657 SRR341660 SRR341661 SRR341663 SRR341664 
# 100       100       100       100       100       100       100       100 
# SRR341665 SRR341669 SRR341670 SRR341673 SRR341674 SRR341675 SRR341676 SRR341681 
# 100       100       100       100       100       100       100       100 
# SRR341684 SRR341687 SRR341693 SRR341713 SRR413575 SRR413576 SRR413578 SRR413579 
# 100       100       100       100       100       100       100       100 
# SRR413580 SRR413581 SRR413582 SRR413584 SRR413585 SRR413587 SRR413592 SRR413593 
# 100       100       100       100       100       100       100       100 
# SRR413594 SRR413597 SRR413598 SRR413599 SRR413600 SRR413601 SRR413603 SRR413605 
# 100       100       100       100       100       100       100       100 
# SRR413606 SRR413607 SRR413608 SRR413610 SRR413613 SRR413614 SRR413615 SRR413616 
# 100       100       100       100       100       100       100       100 
# SRR413617 SRR413618 SRR413619 SRR413620 SRR413621 SRR413623 SRR413625 SRR413626 
# 100       100       100       100       100       100       100       100 
# SRR413634 SRR413637 SRR413642 SRR413652 SRR413660 SRR413661 SRR413670 SRR413688 
# 100       100       100       100       100       100       100       100 
# SRR413758 SRR413768 
# 100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_chn/cpp_analysis_min"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 15270    82
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[81]]
# NULL
# 
# [[82]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_1 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_1 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 SRR341581                 0
# 3 SRR341581                 0
# 4 SRR341581                 0
# 5 SRR341581                 0
# 6 SRR341581                 0
# 7 SRR341581                 0
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# ...
# 
# 
# [1] "completed 80"
# [1] "completed 81"
# [1] "completed 82"
# 
# str(dat)
# 'data.frame':	33618278 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_2" "fxn_2" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5400.094
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 65.8548
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6883
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
# ...
# 
# 
# 
# [[81]]
# NULL
# 
# [[82]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 SRR341581   0.00000000
# 2 cpd02175 SRR341581   0.00000000
# 3 cpd00035 SRR341581   0.13246568
# 4 cpd00117 SRR341581   0.06505299
# 5 cpd00041 SRR341581   0.03059700
# 6 cpd00320 SRR341581   0.00000000
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# ...
# 
# 
# [1] "completed 80"
# [1] "completed 81"
# [1] "completed 82"
# 
# str(dat.cpd.collate)
# 'data.frame':	564406 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1325 0.0651 0.0306 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5400.094
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 65.8548
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund CHN-T2D - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (minimum library size)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-chn-EVEN-sequences/cpp_analysis_min/dat.cpd.collate-all-samps-cpp3d--t2d-chn-rarefied-min-pawsey.rds")

str(dat.cpd.collate)
# 'data.frame':	564406 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1325 0.0651 0.0306 ...

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.000000  0.000000  0.000066  0.009568  0.000959 12.805791 

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf  -4.180    -Inf  -3.018   1.107 


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
# -7.678  -7.678  -4.180  -4.936  -3.018   1.107 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15270 taxa and 82 samples ]
# sample_data() Sample Data:       [ 82 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 15270 taxa by 4 taxonomic ranks ]

head(phy@sam_data)
# Sample Data:        [6 samples by 29 sample variables]:
#   Run actual_read_length..run. Age Assay.Type AvgSpotLen  BioProject    BioSample      Bytes center_name..exp. Center.Name Consent DATASTORE.filetype DATASTORE.provider
# SRR341581 SRR341581                      148  59        WGS        148 PRJNA422434 SAMN00715131 1170792822               BGI         BGI  public          fastq,sra                 s3
# SRR341585 SRR341585                      148  60        WGS        148 PRJNA422434 SAMN00715135 1200370883               BGI         BGI  public          sra,fastq                 s3
# SRR341586 SRR341586                      148  71        WGS        148 PRJNA422434 SAMN00715136 1350066985               BGI         BGI  public          fastq,sra                 s3
# SRR341587 SRR341587                      148  62        WGS        148 PRJNA422434 SAMN00715137 1479828698               BGI         BGI  public          fastq,sra                 s3
# SRR341588 SRR341588                      148  54        WGS        148 PRJNA422434 SAMN00715138 1370794218               BGI         BGI  public          sra,fastq                 s3
# SRR341589 SRR341589                      148  51        WGS        148 PRJNA422434 SAMN00715139 1010953818               BGI         BGI  public          sra,fastq                 s3
# DATASTORE.region Experiment gender                  Instrument     Library.Name LibraryLayout LibrarySelection LibrarySource NATION             Organism Platform
# SRR341581     s3.us-east-1  SRX095662 female Illumina Genome Analyzer II HGMlijMCFDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341585     s3.us-east-1  SRX095666 female Illumina Genome Analyzer II HGMlijMDGDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341586     s3.us-east-1  SRX095667 female Illumina Genome Analyzer II HGMlijMAODFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341587     s3.us-east-1  SRX095668 female Illumina Genome Analyzer II HGMlijMAPDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341588     s3.us-east-1  SRX095669 female Illumina Genome Analyzer II HGMlijMAMDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# SRR341589     s3.us-east-1  SRX095670 female Illumina Genome Analyzer II HGMlijMDSDFAAPEI        PAIRED           RANDOM   METAGENOMIC  China human gut metagenome ILLUMINA
# ReleaseDate   run..run. Sample.Name SRA.Study      Diagnosis
# SRR341581 2012-09-05T00:00:00Z FC615J5AAXX  bgi-DLF001 SRP008047 T2D metformin-
#   SRR341585 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF005 SRP008047 T2D metformin-
#   SRR341586 2012-09-05T00:00:00Z FC61B81AAXX  bgi-DLF006 SRP008047 T2D metformin-
#   SRR341587 2012-09-05T00:00:00Z FC61B81AAXX  bgi-DLF007 SRP008047 T2D metformin-
#   SRR341588 2012-09-05T00:00:00Z FC61B81AAXX  bgi-DLF008 SRP008047 T2D metformin-
#   SRR341589 2012-09-05T00:00:00Z FC61B1KAAXX  bgi-DLF010 SRP008047 T2D metformin-

samp <- as(phy@sam_data,"data.frame")
unique(samp$Diagnosis)
# "T2D metformin-" "ND CTRL"   
samp$group_new <- factor(samp$Diagnosis, 
                         levels = c("T2D metformin-", "ND CTRL"),
                         labels = c("T2D met-", "Normal"),
                         ordered = TRUE )

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( samp$Run )) {
  #i<-1
  this_samp <- samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  dat.cpd.collate$group_label[sel] <- as.character( samp$group_new[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"  
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-min-qty82-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-min-qty82-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	564406 obs. of  5 variables:
# $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0.1325 0.0651 0.0306 ...
# $ log10_abun  : num  -7.678 -7.678 -0.878 -1.187 -1.514 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

length( unique(dat.cpd.collate$cpd_id) ) # 6883
6883*82 # 512387


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 SRR341581   0.00000000 -7.6775653    T2D met-
#   2 cpd02175 SRR341581   0.00000000 -7.6775653    T2D met-
#   3 cpd00035 SRR341581   0.13246568 -0.8778966    T2D met-
#   4 cpd00117 SRR341581   0.06505299 -1.1867327    T2D met-
#   5 cpd00041 SRR341581   0.03059700 -1.5143211    T2D met-
#   6 cpd00320 SRR341581   0.00000000 -7.6775653    T2D met-

dim(data_in) # 564406      5

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

mean(no_compounds) # 4376.024
sd(no_compounds) #  1911.967

mean(sample_sum_relabun) # 65.8548
sd(sample_sum_relabun) # 4.250323

length(unique(data_in$cpd_id)) # 6883

#-------------------------

# 2 of 3 P < 0.05
#### Forslund T2D-CHN - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (minimum library size)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-CHN-T2D-selected-over50s-Host-removal-qty82-EVEN-seqs-min-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-min-qty82-v8e.rds")
str(df) # 'data.frame':	564406 obs. of  5 variables:


## T2D-CHN - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2952 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -7.68 -7.68 -7.68 -7.68 -7.68 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-CHN Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 320, p-value = 4.558e-06
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,8))


p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  ylim( min(res$sum_rel_abun), 0.015 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(a)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-CHN-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-CHN - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-min-qty82-v8e.rds")
str(df) # 'data.frame':	564406 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	410 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd03198" "cpd00224" "cpd00382" "cpd00082" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341581" ...
# $ cpd_rel_abun: num  0.0541 0.5706 0.0446 0.1836 0.0446 ...
# $ log10_abun  : num  -1.267 -0.244 -1.351 -0.736 -1.351 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-CHN Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 1061, p-value = 0.003463
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.58 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(b)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-CHN-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-CHN - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d-ExtraData-Forslund-CHN-T2D-over50s-Hostremoval-EVEN-seqs-min-qty82-v8e.rds")
str(df) # 564406 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	246 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "SRR341581" "SRR341581" "SRR341581" "SRR341585" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -7.68 -7.68 -7.68 -7.68 -7.68 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-CHN Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	82 obs. of  5 variables:
#   $ sample      : chr  "SRR341581" "SRR341585" "SRR341586" "SRR341587" ...
# $ sum_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-CHN Rarefied (Min)" "T2D-CHN Rarefied (Min)" "T2D-CHN Rarefied (Min)" "T2D-CHN Rarefied (Min)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -3.751  -2.782 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.740  -5.740  -5.740  -4.790  -3.751  -2.782 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 643, p-value = 0.07892
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.2 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(c)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-CHN-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


#-------------------------



##########################
########################## RAREFIED SUBSETS ...
########################## T2D-SWE P20
##########################

#### T2D Swedish (SWE) cohort - RERUN subset with even sequences

#### Forslund T2D-SWE - w/ Host-removal - only retain samples with at least >= 20th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-swe.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-swe.rds")

hist(non_host_reads);summary(non_host_reads)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

# only retain samples with at least 1st quartile (>= 25th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.20)
# 20% 
# 5168425 

length(non_host_reads) # 76

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.20)) # 61

keep_t2d_swe_list_20th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_swe_list_20th])
# ERR260231 ERR260139 ERR260140 ERR260206 ERR260234 ERR260203 ERR260244 ERR260255 ERR260169 ERR260205 ERR260227 ERR260147 ERR260207 ERR260153 ERR260241 ERR260210 ERR260258 
# 5168425   5248535   5378909   5548404   5580786   5771403   5982025   6008750   6136558   6410766   6571261   6729275   6768294   6855675   6861705   6950877   7063454 
# ERR260253 ERR260243 ERR260251 ERR260246 ERR260209 ERR260204 ERR260199 ERR260226 ERR260144 ERR260151 ERR260252 ERR260193 ERR260256 ERR260224 ERR260250 ERR260170 ERR260163 
# 7072818   7258793   7403441   7544579   7685166   7771122   7870633   7875104   7906030   7922852   7930550   7966539   7989210   8300315   9039857   9067799   9827775 
# ERR260180 ERR260189 ERR260260 ERR260175 ERR260259 ERR260267 ERR260266 ERR260166 ERR260186 ERR260190 ERR260201 ERR260265 ERR260263 ERR260167 ERR260185 ERR260188 ERR260161 
# 10931504  11522853  11760693  12135382  12156410  12245235  12348901  12817100  13023348  13366593  13576920  13610954  13876773  13934625  13963914  14096899  14235824 
# ERR260165 ERR260198 ERR260264 ERR275252 ERR260179 ERR260181 ERR260174 ERR260162 ERR260171 ERR260173 
# 14392899  14450423  14985725  15782202  16065742  16716732  17868092  19753575  21965572  22466068 

writeLines(keep_t2d_swe_list_20th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/keep_t2d_swe_list_20th.txt")

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 61

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260153"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260151" 
# ERR260153 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_20th/superfocus_out_ERR260153" 

sampid <- keep_t2d_swe_list_20th

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 61


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 454630      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2                           subsys_L3                                                                                                                                fxn percent_abun
# 2 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase                                                                                                      Alanine_racemase_(EC_5.1.1.1) 0.0045927064
# 3 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase                                                                                             Diaminopimelate_epimerase_(EC_5.1.1.7) 0.0021373749
# 4 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase                                                                                                    Glutamate_racemase_(EC_5.1.1.3) 0.0014131404
# 5 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase                           UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0001854747
# 6 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0001854747
# 7 ERR260250 Amino Acids and Derivatives         - Creatine and Creatinine Degradation                                                                                            Creatinine_amidohydrolase_(EC_3.5.2.10) 0.0190773959

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) # 16655    62

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 16655
length(unique(full_fxn_names)) # 16655

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260153"    "ERR260161"    "ERR260162"    "ERR260163"    "ERR260165"    "ERR260166"   
# [13] "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"    "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"    "ERR260180"    "ERR260181"    "ERR260185"    "ERR260186"   
# [25] "ERR260188"    "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"    "ERR260203"    "ERR260204"    "ERR260205"    "ERR260206"    "ERR260207"   
# [37] "ERR260209"    "ERR260210"    "ERR260224"    "ERR260226"    "ERR260227"    "ERR260231"    "ERR260234"    "ERR260241"    "ERR260243"    "ERR260244"    "ERR260246"    "ERR260250"   
# [49] "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"    "ERR260256"    "ERR260258"    "ERR260259"    "ERR260260"    "ERR260263"    "ERR260264"    "ERR260265"    "ERR260266"   
# [61] "ERR260267"    "ERR275252"  

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256"

length(sampid) # 61

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 16655     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 180
length(unique(tax.fxn$subsys_L3)) # 1043
length(unique(tax.fxn$fxn)) # 8918


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 16655     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 16655    61

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16655 taxa and 61 samples ]
# tax_table()   Taxonomy Table:    [ 16655 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260153" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166" "ERR260167" "ERR260169" "ERR260170" "ERR260171" "ERR260173"
# [17] "ERR260174" "ERR260175" "ERR260179" "ERR260180" "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190" "ERR260193" "ERR260198" "ERR260199" "ERR260201" "ERR260203" "ERR260204"
# [33] "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260224" "ERR260226" "ERR260227" "ERR260231" "ERR260234" "ERR260241" "ERR260243" "ERR260244" "ERR260246" "ERR260250" "ERR260251"
# [49] "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260" "ERR260263" "ERR260264" "ERR260265" "ERR260266" "ERR260267" "ERR275252"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse subset of previous fxn phyloseq object
temp <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)
#                Sample Country.subset         Status      Bases       Run group_label      age non_host_reads fxn_sum_counts
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139    T2D met- 70.25205        5248535         198299
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140    T2D met- 70.15342        5378909         255775
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144    T2D met- 69.57534        7906030         268694
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147      Normal 71.39452        6729275         374375
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151    T2D met- 71.56712        7922852         200815
# ERR260153 NG-5636_381            SWE        ND CTRL 2811341262 ERR260153      Normal 70.42466        6855675         283145

# remove columns: 'Bases', non_host_reads, fxn_sum_counts, as not applicable to this version based on rarefied sequences
dim(df.samp) #  61  9
sel <- which(names(df.samp) %in% c("Bases","non_host_reads","fxn_sum_counts"))
df.samp <- df.samp[ ,-sel]
head(df.samp)

# reorder to align with phy object
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16655 taxa and 61 samples ]
# sample_data() Sample Data:       [ 61 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16655 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                                               
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                                         
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                                         
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                                   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"                                                   
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"                                                      
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)"

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

table(phy@sam_data$group_label)
# T2D met-   Normal 
# 28       33 

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 61
NonZeroFxns

mean(NonZeroFxns) # 7452.951
sd(NonZeroFxns) # 1216.816


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 20th percentile)
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
# # For study - Forslund et al T2D-SWE rarefied sequences - 20th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th/working"
# 
# message("\nthis_study <- '-t2d-swe-rarefied-20th-pawsey'")
# this_study <- "-t2d-swe-rarefied-20th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS")
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

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 20th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th/working'
# 
# this_study <- '-t2d-swe-rarefied-20th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS')
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
# [1] 16655     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[16653]]
# NULL
# 
# [[16654]]
# NULL
# 
# [[16655]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 16655
# [1] "added df 1 of 16655"
# [1] "added df 2 of 16655"
# [1] "added df 3 of 16655"
# ...
# 
# 
# [1] "added df 16653 of 16655"
# [1] "added df 16654 of 16655"
# [1] "added df 16655 of 16655"
# 
# str(df.out)
# 'data.frame':	467410 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_1" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id        : chr  NA "rxn03060" "rxn03060" "rxn03060" ...
# $ cpd_id        : chr  NA "cpd01501" "cpd00001" "cpd02597" ...
# $ cpd_name      : chr  NA "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" ...
# $ cpd_form      : chr  NA "C7H7O7" "H2O" "C7H5O6" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                                   f__in
# 1           <NA> NA                                                    <NA>
#   2          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 467409      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9317
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16655 taxa and 61 samples ]
# sample_data() Sample Data:       [ 61 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16655 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 55.94116
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 9315"
# [1] "completed 9316"
# [1] "completed 9317"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9317
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260153 ERR260161 ERR260162 
# 100       100       100       100       100       100       100       100 
# ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 ERR260171 ERR260173 
# 100       100       100       100       100       100       100       100 
# ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 ERR260186 ERR260188 
# 100       100       100       100       100       100       100       100 
# ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 ERR260203 ERR260204 
# 100       100       100       100       100       100       100       100 
# ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260224 ERR260226 ERR260227 
# 100       100       100       100       100       100       100       100 
# ERR260231 ERR260234 ERR260241 ERR260243 ERR260244 ERR260246 ERR260250 ERR260251 
# 100       100       100       100       100       100       100       100 
# ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 ERR260259 ERR260260 ERR260263 
# 100       100       100       100       100       100       100       100 
# ERR260264 ERR260265 ERR260266 ERR260267 ERR275252 
# 100       100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_20th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 16655    61
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# 
# [[59]]
# NULL
# 
# [[60]]
# NULL
# 
# [[61]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                                   f__in
# 2          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 7          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# 7 rxn17391 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm    sample cpd_rel_abun_norm
# 2              1          0.05555556 ERR260139                 0
# 3              1          0.05555556 ERR260139                 0
# 4              1          0.05555556 ERR260139                 0
# 5              1          0.05555556 ERR260139                 0
# 6              1          0.05555556 ERR260139                 0
# 7              1          0.05555556 ERR260139                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 59"
# [1] "completed 60"
# [1] "completed 61"
# 
# str(dat)
# 'data.frame':	28511949 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_1" "fxn_1" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn03060" "rxn03060" "rxn03060" "rxn17391" ...
# $ cpd_id             : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ cpd_name           : chr  "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.0556 0.0556 0.0556 0.0556 0.0556 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4313.498
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.71308
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6913
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
# ...
# 
# 
# 
# [[60]]
# NULL
# 
# [[61]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd01501 ERR260139 0.0000000000
# 2 cpd00001 ERR260139 4.9847302932
# 3 cpd02597 ERR260139 0.0001912333
# 4 cpd24620 ERR260139 0.0000000000
# 5 cpd00851 ERR260139 0.0012344330
# 6 cpd02175 ERR260139 0.0000000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 59"
# [1] "completed 60"
# [1] "completed 61"
# 
# str(dat.cpd.collate)
# 'data.frame':	421693 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.98473 0.000191 0 0.001234 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4313.498
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.71308
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 20th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/cpp_analysis_20th/dat.cpd.collate-all-samps-cpp3d--t2d-swe-rarefied-20th-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000167 0.010229 0.001629 7.156753 

length(unique(dat.cpd.collate$cpd_id)) # 6913
length(unique(dat.cpd.collate$sample)) # 61
str(dat.cpd.collate)
# 'data.frame':	421693 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.98473 0.000191 0 0.001234 ...
6913*61 # 421693

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf -7.1937 -3.7767    -Inf -2.7881  0.8547 

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
# -8.1998 -7.1937 -3.7767 -4.5673 -2.7881  0.8547 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-20th-qty61-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-20th-qty61-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	421693 obs. of  5 variables:
# $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.98473 0.000191 0 0.001234 ...
# $ log10_abun  : num  -8.2 0.698 -3.718 -8.2 -2.909 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd01501 ERR260139 0.0000000000 -8.1998495    T2D met-
# 2 cpd00001 ERR260139 4.9847302932  0.6976417    T2D met-
# 3 cpd02597 ERR260139 0.0001912333 -3.7184364    T2D met-
# 4 cpd24620 ERR260139 0.0000000000 -8.1998495    T2D met-
# 5 cpd00851 ERR260139 0.0012344330 -2.9085325    T2D met-
# 6 cpd02175 ERR260139 0.0000000000 -8.1998495    T2D met-

dim(data_in) #  421693      5

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

mean(no_compounds) # 5191.197
sd(no_compounds) #  497.0133

mean(sample_sum_relabun) # 70.71308
sd(sample_sum_relabun) # 2.912811

length(unique(data_in$cpd_id)) # 6913

#-------------------------


# 2 out of 3 - p < 0.05
#### Forslund T2D-SWE - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 20th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty61-Hostremoval-EVEN-seqs-20th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-20th-qty61-v8e.rds")
str(df) # 'data.frame':	421693 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2196 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  4.38e-06 3.69e-06 3.69e-06 3.69e-06 3.69e-06 ...
# $ log10_abun  : num  -5.36 -5.43 -5.43 -5.43 -5.43 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 28
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 33

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 381, p-value = 0.123
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.0075 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-20th-qty61-v8e.rds")
str(df) # 'data.frame':	421693 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	305 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.1313 0.0798 0.0804 0.2163 0.104 ...
# $ log10_abun  : num  -0.882 -1.098 -1.095 -0.665 -0.983 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 611, p-value = 0.01545
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.56 )+
  expand_limits(y = 1.05*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-20th-qty61-v8e.rds")
str(df) # 'data.frame':	421693 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	183 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -8.2 -8.2 -8.2 -8.2 -8.2 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE Rarefied (P20)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	61 obs. of  5 variables:
#   $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.002783 0.000209 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 2 1 2 1 1 2 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE Rarefied (P20)" "T2D-SWE Rarefied (P20)" "T2D-SWE Rarefied (P20)" "T2D-SWE Rarefied (P20)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -4.561  -2.555 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.387  -6.387  -6.387  -5.402  -4.561  -2.555 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 320.5, p-value = 0.01369
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-20thperc-even-sequences-T2D-SWE-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------



##########################
########################## T2D-SWE P15
##########################
##########################

#### T2D Swedish (SWE) cohort - RERUN subset with even sequences

#### Forslund T2D-SWE - w/ Host-removal - only retain samples with at least >= 15th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-swe.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-swe.rds")

hist(non_host_reads);summary(non_host_reads)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

# only retain samples with at least 1st quartile (>= 15th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.15)
# 15% 
# 5018997 

length(non_host_reads) # 76

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.15)) # 64

keep_t2d_swe_list_15th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_swe_list_15th])
# ERR260230 ERR260225 ERR260159 ERR260231 ERR260139 ERR260140 ERR260206 ERR260234 ERR260203 ERR260244 ERR260255 ERR260169 ERR260205 ERR260227 ERR260147 ERR260207 ERR260153 
# 5038200   5115407   5116818   5168425   5248535   5378909   5548404   5580786   5771403   5982025   6008750   6136558   6410766   6571261   6729275   6768294   6855675 
# ERR260241 ERR260210 ERR260258 ERR260253 ERR260243 ERR260251 ERR260246 ERR260209 ERR260204 ERR260199 ERR260226 ERR260144 ERR260151 ERR260252 ERR260193 ERR260256 ERR260224 
# 6861705   6950877   7063454   7072818   7258793   7403441   7544579   7685166   7771122   7870633   7875104   7906030   7922852   7930550   7966539   7989210   8300315 
# ERR260250 ERR260170 ERR260163 ERR260180 ERR260189 ERR260260 ERR260175 ERR260259 ERR260267 ERR260266 ERR260166 ERR260186 ERR260190 ERR260201 ERR260265 ERR260263 ERR260167 
# 9039857   9067799   9827775  10931504  11522853  11760693  12135382  12156410  12245235  12348901  12817100  13023348  13366593  13576920  13610954  13876773  13934625 
# ERR260185 ERR260188 ERR260161 ERR260165 ERR260198 ERR260264 ERR275252 ERR260179 ERR260181 ERR260174 ERR260162 ERR260171 ERR260173 
# 13963914  14096899  14235824  14392899  14450423  14985725  15782202  16065742  16716732  17868092  19753575  21965572  22466068 

writeLines(keep_t2d_swe_list_15th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/keep_t2d_swe_list_15th.txt")

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 64

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260153"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260151" 
# ERR260153 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_15th/superfocus_out_ERR260153" 

sampid <- keep_t2d_swe_list_15th

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 64


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 476171      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2                           subsys_L3
# 2 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 3 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 4 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 5 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 6 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 7 ERR260250 Amino Acids and Derivatives         - Creatine and Creatinine Degradation
# fxn percent_abun
# 2                                                                                                      Alanine_racemase_(EC_5.1.1.1) 0.0047288031
# 3                                                                                             Diaminopimelate_epimerase_(EC_5.1.1.7) 0.0022007122
# 4                                                                                                    Glutamate_racemase_(EC_5.1.1.3) 0.0014550164
# 5                           UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0001909709
# 6 UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0001909709
# 7                                                                                            Creatinine_amidohydrolase_(EC_3.5.2.10) 0.0190970896

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) #  16664    65

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 16664
length(unique(full_fxn_names)) # 16664

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260153"    "ERR260159"    "ERR260161"    "ERR260162"    "ERR260163"    "ERR260165"   
# [13] "ERR260166"    "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"    "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"    "ERR260180"    "ERR260181"    "ERR260185"   
# [25] "ERR260186"    "ERR260188"    "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"    "ERR260203"    "ERR260204"    "ERR260205"    "ERR260206"   
# [37] "ERR260207"    "ERR260209"    "ERR260210"    "ERR260224"    "ERR260225"    "ERR260226"    "ERR260227"    "ERR260230"    "ERR260231"    "ERR260234"    "ERR260241"    "ERR260243"   
# [49] "ERR260244"    "ERR260246"    "ERR260250"    "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"    "ERR260256"    "ERR260258"    "ERR260259"    "ERR260260"    "ERR260263"   
# [61] "ERR260264"    "ERR260265"    "ERR260266"    "ERR260267"    "ERR275252"   

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256"

length(sampid) # 64

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 16664     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 181
length(unique(tax.fxn$subsys_L3)) # 1047
length(unique(tax.fxn$fxn)) # 8932


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 16664     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 16664    64

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16664 taxa and 64 samples ]
# tax_table()   Taxonomy Table:    [ 16664 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166" "ERR260167" "ERR260169" "ERR260170" "ERR260171"
# [17] "ERR260173" "ERR260174" "ERR260175" "ERR260179" "ERR260180" "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190" "ERR260193" "ERR260198" "ERR260199" "ERR260201" "ERR260203"
# [33] "ERR260204" "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260224" "ERR260225" "ERR260226" "ERR260227" "ERR260230" "ERR260231" "ERR260234" "ERR260241" "ERR260243" "ERR260244"
# [49] "ERR260246" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260" "ERR260263" "ERR260264" "ERR260265" "ERR260266" "ERR260267" "ERR275252"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse subset of previous fxn phyloseq object
temp <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)
#                Sample Country.subset         Status      Bases       Run group_label      age non_host_reads fxn_sum_counts
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139    T2D met- 70.25205        5248535         198299
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140    T2D met- 70.15342        5378909         255775
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144    T2D met- 69.57534        7906030         268694
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147      Normal 71.39452        6729275         374375
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151    T2D met- 71.56712        7922852         200815
# ERR260153 NG-5636_381            SWE        ND CTRL 2811341262 ERR260153      Normal 70.42466        6855675         283145

# remove columns: 'Bases', non_host_reads, fxn_sum_counts, as not applicable to this version based on rarefied sequences
dim(df.samp) #  64  9
sel <- which(names(df.samp) %in% c("Bases","non_host_reads","fxn_sum_counts"))
df.samp <- df.samp[ ,-sel]
head(df.samp)

# reorder to align with phy object
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16664 taxa and 64 samples ]
# sample_data() Sample Data:       [ 64 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16664 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                                               
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                                         
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                                         
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                                   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"                                                   
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"                                                      
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)"

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

table(phy@sam_data$group_label)
# T2D met-   Normal 
# 29       35 

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 64
NonZeroFxns

mean(NonZeroFxns) # 7440.172
sd(NonZeroFxns) # 1226.594


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 15th percentile)
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
# # For study - Forslund et al T2D-SWE rarefied sequences - 15th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th/working"
# 
# message("\nthis_study <- '-t2d-swe-rarefied-15th-pawsey'")
# this_study <- "-t2d-swe-rarefied-15th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS")
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

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 15th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th/working'
# 
# this_study <- '-t2d-swe-rarefied-15th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS')
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
# [1] 16664     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# [[16663]]
# NULL
# 
# [[16664]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 16664
# [1] "added df 1 of 16664"
# [1] "added df 2 of 16664"
# [1] "added df 3 of 16664"
# ...
# 
# 
# [1] "added df 16662 of 16664"
# [1] "added df 16663 of 16664"
# [1] "added df 16664 of 16664"
# 
# str(df.out)
# 'data.frame':	459381 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_1" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id        : chr  NA "rxn03060" "rxn03060" "rxn03060" ...
# $ cpd_id        : chr  NA "cpd01501" "cpd00001" "cpd02597" ...
# $ cpd_name      : chr  NA "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" ...
# $ cpd_form      : chr  NA "C7H7O7" "H2O" "C7H5O6" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                                   f__in
# 1           <NA> NA                                                    <NA>
#   2          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 459380      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9318
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16664 taxa and 64 samples ]
# sample_data() Sample Data:       [ 64 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16664 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 55.91695
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 9316"
# [1] "completed 9317"
# [1] "completed 9318"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9318
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260153 ERR260159 ERR260161 
# 100       100       100       100       100       100       100       100 
# ERR260162 ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 ERR260171 
# 100       100       100       100       100       100       100       100 
# ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 ERR260186 
# 100       100       100       100       100       100       100       100 
# ERR260188 ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 ERR260203 
# 100       100       100       100       100       100       100       100 
# ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260224 ERR260225 
# 100       100       100       100       100       100       100       100 
# ERR260226 ERR260227 ERR260230 ERR260231 ERR260234 ERR260241 ERR260243 ERR260244 
# 100       100       100       100       100       100       100       100 
# ERR260246 ERR260250 ERR260251 ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 
# 100       100       100       100       100       100       100       100 
# ERR260259 ERR260260 ERR260263 ERR260264 ERR260265 ERR260266 ERR260267 ERR275252 
# 100       100       100       100       100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_15th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 16664    64
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[63]]
# NULL
# 
# [[64]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                                   f__in
# 2          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 7          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# 7 rxn17391 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm    sample cpd_rel_abun_norm
# 2              1          0.05555556 ERR260139                 0
# 3              1          0.05555556 ERR260139                 0
# 4              1          0.05555556 ERR260139                 0
# 5              1          0.05555556 ERR260139                 0
# 6              1          0.05555556 ERR260139                 0
# 7              1          0.05555556 ERR260139                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 62"
# [1] "completed 63"
# [1] "completed 64"
# 
# str(dat)
# 'data.frame':	29400320 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_1" "fxn_1" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn03060" "rxn03060" "rxn03060" "rxn17391" ...
# $ cpd_id             : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ cpd_name           : chr  "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.0556 0.0556 0.0556 0.0556 0.0556 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4517.961
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.59315
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6930
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
# ...
# 
# 
# 
# [[63]]
# NULL
# 
# [[64]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd01501 ERR260139 0.0000000000
# 2 cpd00001 ERR260139 4.9849481841
# 3 cpd02597 ERR260139 0.0001968933
# 4 cpd24620 ERR260139 0.0000000000
# 5 cpd00851 ERR260139 0.0011956843
# 6 cpd02175 ERR260139 0.0000000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# [1] "completed 63"
# [1] "completed 64"
# 
# str(dat.cpd.collate)
# 'data.frame':	443520 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.984948 0.000197 0 0.001196 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4517.961
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.59315
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 15th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/cpp_analysis_15th/dat.cpd.collate-all-samps-cpp3d--t2d-swe-rarefied-15th-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000167 0.010187 0.001617 7.157112 

length(unique(dat.cpd.collate$cpd_id)) # 6930
length(unique(dat.cpd.collate$sample)) # 64
str(dat.cpd.collate)
# 'data.frame':	443520 obs. of  3 variables:
# $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.984948 0.000197 0 0.001196 ...
6930*64 # 443520

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -3.7776    -Inf -2.7912  0.8547 

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
# -8.1871 -8.1871 -3.7776 -4.5739 -2.7912  0.8547 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-15th-qty64-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-15th-qty64-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	443520 obs. of  5 variables:
# $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.984948 0.000197 0 0.001196 ...
# $ log10_abun  : num  -8.187 0.698 -3.706 -8.187 -2.922 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd01501 ERR260139 0.0000000000 -8.1871438    T2D met-
# 2 cpd00001 ERR260139 4.9849481841  0.6976606    T2D met-
# 3 cpd02597 ERR260139 0.0001968933 -3.7057691    T2D met-
# 4 cpd24620 ERR260139 0.0000000000 -8.1871438    T2D met-
# 5 cpd00851 ERR260139 0.0011956843 -2.9223835    T2D met-
# 6 cpd02175 ERR260139 0.0000000000 -8.1871438    T2D met-

dim(data_in) # 443520      5

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

mean(no_compounds) # 5176.812
sd(no_compounds) #  504.5694

mean(sample_sum_relabun) # 70.59315
sd(sample_sum_relabun) # 2.92125

length(unique(data_in$cpd_id)) # 6930

#-------------------------


# 2 out of 3 P < 0.05
#### Forslund T2D-SWE - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 15th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty64-Hostremoval-EVEN-seqs-15th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-15th-qty64-v8e.rds")
str(df) # 'data.frame':	443520 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2304 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  4.51e-06 3.79e-06 3.79e-06 3.79e-06 3.79e-06 ...
# $ log10_abun  : num  -5.35 -5.42 -5.42 -5.42 -5.42 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 29
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 35

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 413, p-value = 0.1032
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.0075 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-15th-qty64-v8e.rds")
str(df) # 'data.frame':	443520 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	320 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00224" "cpd03198" "cpd00382" "cpd00082" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.1314 0.0799 0.0805 0.2165 0.1045 ...
# $ log10_abun  : num  -0.881 -1.097 -1.094 -0.665 -0.981 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 681, p-value = 0.00945
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.565 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-15th-qty64-v8e.rds")
str(df) # 443520 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	192 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -8.19 -8.19 -8.19 -8.19 -8.19 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE Rarefied (P15)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	64 obs. of  5 variables:
# $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.002805 0.000167 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 2 1 2 1 1 1 2 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE Rarefied (P15)" "T2D-SWE Rarefied (P15)" "T2D-SWE Rarefied (P15)" "T2D-SWE Rarefied (P15)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -4.477  -2.552 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.374  -6.374  -6.374  -5.387  -4.477  -2.552 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 344, p-value = 0.008559
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-15thperc-even-sequences-T2D-SWE-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------




##########################
########################## T2D-SWE P10
##########################
##########################


#### T2D Swedish (SWE) cohort - RERUN subset with even sequences

#### Forslund T2D-SWE - w/ Host-removal - only retain samples with at least >= 10th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-swe.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-swe.rds")

hist(non_host_reads);summary(non_host_reads)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

# only retain samples with at least 1st quartile (>= 10th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.10)
# 10% 
# 4330969 

length(non_host_reads) # 76

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.10)) # 68

keep_t2d_swe_list_10th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_swe_list_10th])
# ERR260215 ERR260152 ERR260242 ERR260217 ERR260230 ERR260225 ERR260159 ERR260231 ERR260139 ERR260140 ERR260206 ERR260234 ERR260203 ERR260244 ERR260255 ERR260169 ERR260205 ERR260227 
# 4469057   4516029   4780520   5012596   5038200   5115407   5116818   5168425   5248535   5378909   5548404   5580786   5771403   5982025   6008750   6136558   6410766   6571261 
# ERR260147 ERR260207 ERR260153 ERR260241 ERR260210 ERR260258 ERR260253 ERR260243 ERR260251 ERR260246 ERR260209 ERR260204 ERR260199 ERR260226 ERR260144 ERR260151 ERR260252 ERR260193 
# 6729275   6768294   6855675   6861705   6950877   7063454   7072818   7258793   7403441   7544579   7685166   7771122   7870633   7875104   7906030   7922852   7930550   7966539 
# ERR260256 ERR260224 ERR260250 ERR260170 ERR260163 ERR260180 ERR260189 ERR260260 ERR260175 ERR260259 ERR260267 ERR260266 ERR260166 ERR260186 ERR260190 ERR260201 ERR260265 ERR260263 
# 7989210   8300315   9039857   9067799   9827775  10931504  11522853  11760693  12135382  12156410  12245235  12348901  12817100  13023348  13366593  13576920  13610954  13876773 
# ERR260167 ERR260185 ERR260188 ERR260161 ERR260165 ERR260198 ERR260264 ERR275252 ERR260179 ERR260181 ERR260174 ERR260162 ERR260171 ERR260173 
# 13934625  13963914  14096899  14235824  14392899  14450423  14985725  15782202  16065742  16716732  17868092  19753575  21965572  22466068 

writeLines(keep_t2d_swe_list_10th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/keep_t2d_swe_list_10th.txt")

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 68

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260152"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260151" 
# ERR260152 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_10th/superfocus_out_ERR260152" 

sampid <- keep_t2d_swe_list_10th

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 68


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 491069      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2                           subsys_L3
# 2 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 3 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 4 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 5 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 6 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 7 ERR260250 Amino Acids and Derivatives         - Creatine and Creatinine Degradation
# fxn percent_abun
# 2                                                                                                      Alanine_racemase_(EC_5.1.1.1) 0.0040002105
# 3                                                                                             Diaminopimelate_epimerase_(EC_5.1.1.7) 0.0024211801
# 4                                                                                                    Glutamate_racemase_(EC_5.1.1.3) 0.0010526870
# 5                           UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0002210643
# 6 UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0002210643
# 7                                                                                            Creatinine_amidohydrolase_(EC_3.5.2.10) 0.0195799779

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) #  16485    69

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 16485
length(unique(full_fxn_names)) # 16485

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Aspartate_racemase_(EC_5.1.1.13)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260152"    "ERR260153"    "ERR260159"    "ERR260161"    "ERR260162"    "ERR260163"   
# [13] "ERR260165"    "ERR260166"    "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"    "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"    "ERR260180"    "ERR260181"   
# [25] "ERR260185"    "ERR260186"    "ERR260188"    "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"    "ERR260203"    "ERR260204"    "ERR260205"   
# [37] "ERR260206"    "ERR260207"    "ERR260209"    "ERR260210"    "ERR260215"    "ERR260217"    "ERR260224"    "ERR260225"    "ERR260226"    "ERR260227"    "ERR260230"    "ERR260231"   
# [49] "ERR260234"    "ERR260241"    "ERR260242"    "ERR260243"    "ERR260244"    "ERR260246"    "ERR260250"    "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"    "ERR260256"   
# [61] "ERR260258"    "ERR260259"    "ERR260260"    "ERR260263"    "ERR260264"    "ERR260265"    "ERR260266"    "ERR260267"    "ERR275252"   

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256"

length(sampid) # 68

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 16485     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 181
length(unique(tax.fxn$subsys_L3)) # 1043
length(unique(tax.fxn$fxn)) # 8863


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 16485     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 16485    68

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16485 taxa and 68 samples ]
# tax_table()   Taxonomy Table:    [ 16485 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166" "ERR260167" "ERR260169" "ERR260170"
# [17] "ERR260171" "ERR260173" "ERR260174" "ERR260175" "ERR260179" "ERR260180" "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190" "ERR260193" "ERR260198" "ERR260199" "ERR260201"
# [33] "ERR260203" "ERR260204" "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260215" "ERR260217" "ERR260224" "ERR260225" "ERR260226" "ERR260227" "ERR260230" "ERR260231" "ERR260234"
# [49] "ERR260241" "ERR260242" "ERR260243" "ERR260244" "ERR260246" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260" "ERR260263" "ERR260264"
# [65] "ERR260265" "ERR260266" "ERR260267" "ERR275252"

### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse subset of previous fxn phyloseq object
temp <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)
#                Sample Country.subset         Status      Bases       Run group_label      age non_host_reads fxn_sum_counts
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139    T2D met- 70.25205        5248535         198299
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140    T2D met- 70.15342        5378909         255775
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144    T2D met- 69.57534        7906030         268694
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147      Normal 71.39452        6729275         374375
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151    T2D met- 71.56712        7922852         200815
# ERR260152 NG-5636_380            SWE T2D metformin- 1813559434 ERR260152    T2D met- 71.24384        4516029         154717

# remove columns: 'Bases', non_host_reads, fxn_sum_counts, as not applicable to this version based on rarefied sequences
dim(df.samp) #  68  9
sel <- which(names(df.samp) %in% c("Bases","non_host_reads","fxn_sum_counts"))
df.samp <- df.samp[ ,-sel]
head(df.samp)

# reorder to align with phy object
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16485 taxa and 68 samples ]
# sample_data() Sample Data:       [ 68 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16485 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                                               
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                                         
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                                   
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"                                                   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"                                                      
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)"
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Aspartate_racemase_(EC_5.1.1.13)"                            

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

table(phy@sam_data$group_label)
# T2D met-   Normal 
# 30       38 

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 68
NonZeroFxns

mean(NonZeroFxns) # 7221.603
sd(NonZeroFxns) # 1186.78


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 10th percentile)
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
# # For study - Forslund et al T2D-SWE rarefied sequences - 10th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th/working"
# 
# message("\nthis_study <- '-t2d-swe-rarefied-10th-pawsey'")
# this_study <- "-t2d-swe-rarefied-10th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS")
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

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 10th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th/working'
# 
# this_study <- '-t2d-swe-rarefied-10th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS')
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
# [1] 16485     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[16483]]
# NULL
# 
# [[16484]]
# NULL
# 
# [[16485]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 16485
# [1] "added df 1 of 16485"
# [1] "added df 2 of 16485"
# [1] "added df 3 of 16485"
# ...
# 
# 
# [1] "added df 16483 of 16485"
# [1] "added df 16484 of 16485"
# [1] "added df 16485 of 16485"
# 
# str(df.out)
# 'data.frame':	454178 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_1" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id        : chr  NA "rxn03060" "rxn03060" "rxn03060" ...
# $ cpd_id        : chr  NA "cpd01501" "cpd00001" "cpd02597" ...
# $ cpd_name      : chr  NA "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" ...
# $ cpd_form      : chr  NA "C7H7O7" "H2O" "C7H5O6" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                                   f__in
# 1           <NA> NA                                                    <NA>
#   2          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 454177      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 9222
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16485 taxa and 68 samples ]
# sample_data() Sample Data:       [ 68 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 16485 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 55.94177
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 9220"
# [1] "completed 9221"
# [1] "completed 9222"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 9222
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260152 ERR260153 ERR260159 
# 100       100       100       100       100       100       100       100 
# ERR260161 ERR260162 ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 
# 100       100       100       100       100       100       100       100 
# ERR260171 ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 
# 100       100       100       100       100       100       100       100 
# ERR260186 ERR260188 ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 
# 100       100       100       100       100       100       100       100 
# ERR260203 ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260215 
# 100       100       100       100       100       100       100       100 
# ERR260217 ERR260224 ERR260225 ERR260226 ERR260227 ERR260230 ERR260231 ERR260234 
# 100       100       100       100       100       100       100       100 
# ERR260241 ERR260242 ERR260243 ERR260244 ERR260246 ERR260250 ERR260251 ERR260252 
# 100       100       100       100       100       100       100       100 
# ERR260253 ERR260255 ERR260256 ERR260258 ERR260259 ERR260260 ERR260263 ERR260264 
# 100       100       100       100       100       100       100       100 
# ERR260265 ERR260266 ERR260267 ERR275252 
# 100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_10th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 16485    68
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[66]]
# NULL
# 
# [[67]]
# NULL
# 
# [[68]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                                   f__in
# 2          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 7          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# 7 rxn17391 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm    sample cpd_rel_abun_norm
# 2              1          0.05555556 ERR260139                 0
# 3              1          0.05555556 ERR260139                 0
# 4              1          0.05555556 ERR260139                 0
# 5              1          0.05555556 ERR260139                 0
# 6              1          0.05555556 ERR260139                 0
# 7              1          0.05555556 ERR260139                 0
# [1] "completed 2"
# [1] "completed 3"
# [1] "completed 4"
# ...
# 
# 
# [1] "completed 66"
# [1] "completed 67"
# [1] "completed 68"
# 
# str(dat)
# 'data.frame':	30884036 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_1" "fxn_1" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn03060" "rxn03060" "rxn03060" "rxn17391" ...
# $ cpd_id             : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ cpd_name           : chr  "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.0556 0.0556 0.0556 0.0556 0.0556 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 4805.322
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.6665
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6927
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
# ...
# 
# 
# 
# [[67]]
# NULL
# 
# [[68]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd01501 ERR260139 0.0000000000
# 2 cpd00001 ERR260139 4.9802357967
# 3 cpd02597 ERR260139 0.0001129524
# 4 cpd24620 ERR260139 0.0000000000
# 5 cpd00035 ERR260139 0.2561419336
# 6 cpd00117 ERR260139 0.2054878111
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 66"
# [1] "completed 67"
# [1] "completed 68"
# 
# str(dat.cpd.collate)
# 'data.frame':	471036 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.980236 0.000113 0 0.256142 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 4805.322
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.6665
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 10th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/cpp_analysis_10th/dat.cpd.collate-all-samps-cpp3d--t2d-swe-rarefied-10th-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000165 0.010202 0.001615 7.160229 

length(unique(dat.cpd.collate$cpd_id)) # 6927
length(unique(dat.cpd.collate$sample)) # 68
str(dat.cpd.collate)
# 'data.frame':	471036 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.980236 0.000113 0 0.256142 ...
6927*68 # 471036

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -3.7819    -Inf -2.7919  0.8549 

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
# -8.1236 -8.1236 -3.7819 -4.5800 -2.7919  0.8549 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-10th-qty68-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-10th-qty68-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	471036 obs. of  5 variables:
# $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 4.980236 0.000113 0 0.256142 ...
# $ log10_abun  : num  -8.124 0.697 -3.947 -8.124 -0.592 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd01501 ERR260139 0.0000000000 -8.1235624    T2D met-
#   2 cpd00001 ERR260139 4.9802357967  0.6972499    T2D met-
#   3 cpd02597 ERR260139 0.0001129524 -3.9471045    T2D met-
#   4 cpd24620 ERR260139 0.0000000000 -8.1235624    T2D met-
#   5 cpd00035 ERR260139 0.2561419336 -0.5915193    T2D met-
#   6 cpd00117 ERR260139 0.2054878111 -0.6872139    T2D met-

dim(data_in) # 471036      5

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

mean(no_compounds) # 5119.412
sd(no_compounds) #  501.2984

mean(sample_sum_relabun) # 70.6665
sd(sample_sum_relabun) # 2.866425

length(unique(data_in$cpd_id)) # 6927

#-------------------------


# 2 of 3 P < 0.05
#### Forslund T2D-SWE - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 10th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty68-Hostremoval-EVEN-seqs-10th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-10th-qty68-v8e.rds")
str(df) # 'data.frame':	471036 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2448 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  4.81e-06 4.40e-06 4.40e-06 4.40e-06 4.40e-06 ...
# $ log10_abun  : num  -5.32 -5.36 -5.36 -5.36 -5.36 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 30
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 38

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 441, p-value = 0.05644
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  ylim( min(res$sum_rel_abun), 0.0075 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-10th-qty68-v8e.rds")
str(df) # 'data.frame':	471036 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	340 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd03198" "cpd00224" "cpd00382" "cpd00082" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.0798 0.1321 0.0804 0.2182 0.1053 ...
# $ log10_abun  : num  -1.098 -0.879 -1.095 -0.661 -0.978 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 761, p-value = 0.008964
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.61 )+
  expand_limits(y = 1.15*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-10th-qty68-v8e.rds")
str(df) # 471036 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	204 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -8.12 -8.12 -8.12 -8.12 -8.12 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE Rarefied (P10)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	68 obs. of  5 variables:
#   $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.002807 0.000134 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 2 1 1 2 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE Rarefied (P10)" "T2D-SWE Rarefied (P10)" "T2D-SWE Rarefied (P10)" "T2D-SWE Rarefied (P10)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -4.479  -2.552 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -6.311  -6.311  -6.311  -5.363  -4.479  -2.552 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 385.5, p-value = 0.006897
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-10thperc-even-sequences-T2D-SWE-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------




##########################
########################## T2D-SWE P5
##########################
##########################


#### T2D Swedish (SWE) cohort - RERUN subset with even sequences

#### Forslund T2D-SWE - w/ Host-removal - only retain samples with at least >= 5th percentile number of sequences
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-swe.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-swe.rds")

hist(non_host_reads);summary(non_host_reads)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

sum(non_host_reads)

# only retain samples with at least 1st quartile (>= 5th percentile) number of sequences

quantile(x = non_host_reads, probs = 0.05)
# 5% 
# 2454297 

length(non_host_reads) # 72

sel <- which(non_host_reads >= quantile(x = non_host_reads, probs = 0.05)) # 72

keep_t2d_swe_list_5th <- names(non_host_reads)[sel]

sort( non_host_reads[keep_t2d_swe_list_5th])
# ERR260216 ERR260223 ERR260218 ERR260221 ERR260215 ERR260152 ERR260242 ERR260217 ERR260230 ERR260225 ERR260159 ERR260231 ERR260139 ERR260140 ERR260206 ERR260234 ERR260203 
# 2630853   3349000   4118993   4192881   4469057   4516029   4780520   5012596   5038200   5115407   5116818   5168425   5248535   5378909   5548404   5580786   5771403 
# ERR260244 ERR260255 ERR260169 ERR260205 ERR260227 ERR260147 ERR260207 ERR260153 ERR260241 ERR260210 ERR260258 ERR260253 ERR260243 ERR260251 ERR260246 ERR260209 ERR260204 
# 5982025   6008750   6136558   6410766   6571261   6729275   6768294   6855675   6861705   6950877   7063454   7072818   7258793   7403441   7544579   7685166   7771122 
# ERR260199 ERR260226 ERR260144 ERR260151 ERR260252 ERR260193 ERR260256 ERR260224 ERR260250 ERR260170 ERR260163 ERR260180 ERR260189 ERR260260 ERR260175 ERR260259 ERR260267 
# 7870633   7875104   7906030   7922852   7930550   7966539   7989210   8300315   9039857   9067799   9827775  10931504  11522853  11760693  12135382  12156410  12245235 
# ERR260266 ERR260166 ERR260186 ERR260190 ERR260201 ERR260265 ERR260263 ERR260167 ERR260185 ERR260188 ERR260161 ERR260165 ERR260198 ERR260264 ERR275252 ERR260179 ERR260181 
# 12348901  12817100  13023348  13366593  13576920  13610954  13876773  13934625  13963914  14096899  14235824  14392899  14450423  14985725  15782202  16065742  16716732 
# ERR260174 ERR260162 ERR260171 ERR260173 
# 17868092  19753575  21965572  22466068 

writeLines(keep_t2d_swe_list_5th, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/keep_t2d_swe_list_5th.txt")

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 72

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260152"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260151" 
# ERR260152 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_5th/superfocus_out_ERR260152" 

sampid <- keep_t2d_swe_list_5th

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 72


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 474016      6
head(sfx.long)
# sampleID                   subsys_L1 subsys_L2                           subsys_L3
# 2 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 3 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 4 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 5 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 6 ERR260250 Amino Acids and Derivatives         -                 Amino acid racemase
# 7 ERR260250 Amino Acids and Derivatives         - Creatine and Creatinine Degradation
# fxn percent_abun
# 2                                                                                                      Alanine_racemase_(EC_5.1.1.1) 0.0033380808
# 3                                                                                             Diaminopimelate_epimerase_(EC_5.1.1.7) 0.0021326627
# 4                                                                                                    Glutamate_racemase_(EC_5.1.1.3) 0.0007417957
# 5                           UDP-N-acetylmuramoyl-tripeptide--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0003894428
# 6 UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine_ligase_(EC_6.3.2.10)_/_Alanine_racemase_(EC_5.1.1.1) 0.0001112694
# 7                                                                                            Creatinine_amidohydrolase_(EC_3.5.2.10) 0.0178030977

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) #  15859    73

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 15859
length(unique(full_fxn_names)) # 15859

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260152"    "ERR260153"    "ERR260159"    "ERR260161"    "ERR260162"   
# [12] "ERR260163"    "ERR260165"    "ERR260166"    "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"    "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"   
# [23] "ERR260180"    "ERR260181"    "ERR260185"    "ERR260186"    "ERR260188"    "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"   
# [34] "ERR260203"    "ERR260204"    "ERR260205"    "ERR260206"    "ERR260207"    "ERR260209"    "ERR260210"    "ERR260215"    "ERR260216"    "ERR260217"    "ERR260218"   
# [45] "ERR260221"    "ERR260223"    "ERR260224"    "ERR260225"    "ERR260226"    "ERR260227"    "ERR260230"    "ERR260231"    "ERR260234"    "ERR260241"    "ERR260242"   
# [56] "ERR260243"    "ERR260244"    "ERR260246"    "ERR260250"    "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"    "ERR260256"    "ERR260258"    "ERR260259"   
# [67] "ERR260260"    "ERR260263"    "ERR260264"    "ERR260265"    "ERR260266"    "ERR260267"    "ERR275252"   

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256"

length(sampid) # 72

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 15859     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 178
length(unique(tax.fxn$subsys_L3)) # 1017
length(unique(tax.fxn$fxn)) # 8581


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 15859     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 15859    72

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15859 taxa and 72 samples ]
# tax_table()   Taxonomy Table:    [ 15859 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166" "ERR260167"
# [15] "ERR260169" "ERR260170" "ERR260171" "ERR260173" "ERR260174" "ERR260175" "ERR260179" "ERR260180" "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190"
# [29] "ERR260193" "ERR260198" "ERR260199" "ERR260201" "ERR260203" "ERR260204" "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260215" "ERR260216" "ERR260217"
# [43] "ERR260218" "ERR260221" "ERR260223" "ERR260224" "ERR260225" "ERR260226" "ERR260227" "ERR260230" "ERR260231" "ERR260234" "ERR260241" "ERR260242" "ERR260243" "ERR260244"
# [57] "ERR260246" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260" "ERR260263" "ERR260264" "ERR260265" "ERR260266"
# [71] "ERR260267" "ERR275252"


### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse subset of previous fxn phyloseq object
temp <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)
#                Sample Country.subset         Status      Bases       Run group_label      age non_host_reads fxn_sum_counts
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139    T2D met- 70.25205        5248535         198299
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140    T2D met- 70.15342        5378909         255775
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144    T2D met- 69.57534        7906030         268694
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147      Normal 71.39452        6729275         374375
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151    T2D met- 71.56712        7922852         200815
# ERR260152 NG-5636_380            SWE T2D metformin- 1813559434 ERR260152    T2D met- 71.24384        4516029         154717

# remove columns: 'Bases', non_host_reads, fxn_sum_counts, as not applicable to this version based on rarefied sequences
dim(df.samp) #  72  9
sel <- which(names(df.samp) %in% c("Bases","non_host_reads","fxn_sum_counts"))
df.samp <- df.samp[ ,-sel]
head(df.samp)

# reorder to align with phy object
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15859 taxa and 72 samples ]
# sample_data() Sample Data:       [ 72 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 15859 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                                               
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "2-methylcitrate_dehydratase_FeS_dependent_(EC_4.2.1.79)"                                         
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                                         
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                                   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"                                                   
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"                                                      
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)"                       

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

table(phy@sam_data$group_label)
# T2D met-   Normal 
# 30       42 

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 72
NonZeroFxns

mean(NonZeroFxns) # 6583.556
sd(NonZeroFxns) # 1145.311


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (>= 5th percentile)
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
# # For study - Forslund et al T2D-SWE rarefied sequences - 5th percentile
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th/working"
# 
# message("\nthis_study <- '-t2d-swe-rarefied-5th-pawsey'")
# this_study <- "-t2d-swe-rarefied-5th-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS")
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

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (>= 5th percentile)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th/working'
# 
# this_study <- '-t2d-swe-rarefied-5th-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS')
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
# [1] 15859     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[15857]]
# NULL
# 
# [[15858]]
# NULL
# 
# [[15859]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 15859
# [1] "added df 1 of 15859"
# [1] "added df 2 of 15859"
# [1] "added df 3 of 15859"
# ...
# 
# [1] "added df 15856 of 15859"
# [1] "added df 15857 of 15859"
# [1] "added df 15858 of 15859"
# [1] "added df 15859 of 15859"
# 
# str(df.out)
# 'data.frame':	429863 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_1" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id        : chr  NA "rxn03060" "rxn03060" "rxn03060" ...
# $ cpd_id        : chr  NA "cpd01501" "cpd00001" "cpd02597" ...
# $ cpd_name      : chr  NA "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" ...
# $ cpd_form      : chr  NA "C7H7O7" "H2O" "C7H5O6" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                                   f__in
# 1           <NA> NA                                                    <NA>
#   2          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1  1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 1     <NA>     <NA>                                         <NA>     <NA>
#   2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# cpd_molar_prop
# 1             NA
# 2              1
# 3              1
# 4              1
# 5              1
# 6              1
# 
# dim(df.out)
# [1] 429862      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 8898
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 15859 taxa and 72 samples ]
# sample_data() Sample Data:       [ 72 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 15859 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 56.10694
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 8895"
# [1] "completed 8896"
# [1] "completed 8897"
# [1] "completed 8898"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 8898
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260152 ERR260153 ERR260159 
# 100       100       100       100       100       100       100       100 
# ERR260161 ERR260162 ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 
# 100       100       100       100       100       100       100       100 
# ERR260171 ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 
# 100       100       100       100       100       100       100       100 
# ERR260186 ERR260188 ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 
# 100       100       100       100       100       100       100       100 
# ERR260203 ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260215 
# 100       100       100       100       100       100       100       100 
# ERR260216 ERR260217 ERR260218 ERR260221 ERR260223 ERR260224 ERR260225 ERR260226 
# 100       100       100       100       100       100       100       100 
# ERR260227 ERR260230 ERR260231 ERR260234 ERR260241 ERR260242 ERR260243 ERR260244 
# 100       100       100       100       100       100       100       100 
# ERR260246 ERR260250 ERR260251 ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 
# 100       100       100       100       100       100       100       100 
# ERR260259 ERR260260 ERR260263 ERR260264 ERR260265 ERR260266 ERR260267 ERR275252 
# 100       100       100       100       100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_5th"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 15859    72
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
# [[71]]
# NULL
# 
# [[72]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                                   f__in
# 2          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 3          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 4          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 5          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 6          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# 7          fxn_1 1 2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)
# rxn_id   cpd_id                                     cpd_name cpd_form
# 2 rxn03060 cpd01501                              2-Methylcitrate   C7H7O7
# 3 rxn03060 cpd00001                                          H2O      H2O
# 4 rxn03060 cpd02597                        cis-2-Methylaconitate   C7H5O6
# 5 rxn17391 cpd24620 (2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate   C7H7O7
# 6 rxn17391 cpd00001                                          H2O      H2O
# 7 rxn17391 cpd02597                        cis-2-Methylaconitate   C7H5O6
# cpd_molar_prop cpd_molar_prop_norm    sample cpd_rel_abun_norm
# 2              1          0.05555556 ERR260139                 0
# 3              1          0.05555556 ERR260139                 0
# 4              1          0.05555556 ERR260139                 0
# 5              1          0.05555556 ERR260139                 0
# 6              1          0.05555556 ERR260139                 0
# 7              1          0.05555556 ERR260139                 0
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 70"
# [1] "completed 71"
# [1] "completed 72"
# 
# str(dat)
# 'data.frame':	30950064 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_1" "fxn_1" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" "2-methylcitrate dehydratase FeS dependent (EC 4.2.1.79)" ...
# $ rxn_id             : chr  "rxn03060" "rxn03060" "rxn03060" "rxn17391" ...
# $ cpd_id             : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ cpd_name           : chr  "2-Methylcitrate" "H2O" "cis-2-Methylaconitate" "(2S,3S)-2-hydroxybutane-1,2,3-tricarboxylate" ...
# $ cpd_form           : chr  "C7H7O7" "H2O" "C7H5O6" "C7H7O7" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.0556 0.0556 0.0556 0.0556 0.0556 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0 0 0 0 0 0 0 0 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5073.657
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.46745
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6839
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
# ...
# 
# 
# 
# [[71]]
# NULL
# 
# [[72]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd01501 ERR260139 0.0000000000
# 2 cpd00001 ERR260139 5.0015412625
# 3 cpd02597 ERR260139 0.0002028151
# 4 cpd24620 ERR260139 0.0000000000
# 5 cpd00851 ERR260139 0.0011454514
# 6 cpd02175 ERR260139 0.0000000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 70"
# [1] "completed 71"
# [1] "completed 72"
# 
# str(dat.cpd.collate)
# 'data.frame':	492408 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 5.001541 0.000203 0 0.001145 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5073.657
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.46745
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (>= 5th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/cpp_analysis_5th/dat.cpd.collate-all-samps-cpp3d--t2d-swe-rarefied-5th-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000170 0.010304 0.001661 7.137497 

length(unique(dat.cpd.collate$cpd_id)) # 6839
length(unique(dat.cpd.collate$sample)) # 72
str(dat.cpd.collate)
# 'data.frame':	492408 obs. of  3 variables:
# $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 5.001541 0.000203 0 0.001145 ...
6839*72 # 492408

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -3.7686    -Inf -2.7797  0.8535 

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
# -8.1901 -8.1901 -3.7686 -4.6060 -2.7797  0.8535 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-5th-qty72-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-5th-qty72-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	492408 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd01501" "cpd00001" "cpd02597" "cpd24620" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 5.001541 0.000203 0 0.001145 ...
# $ log10_abun  : num  -8.19 0.699 -3.693 -8.19 -2.941 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd01501 ERR260139 0.0000000000 -8.1901366    T2D met-
# 2 cpd00001 ERR260139 5.0015412625  0.6991039    T2D met-
# 3 cpd02597 ERR260139 0.0002028151 -3.6928996    T2D met-
# 4 cpd24620 ERR260139 0.0000000000 -8.1901366    T2D met-
# 5 cpd00851 ERR260139 0.0011454514 -2.9410233    T2D met-
# 6 cpd02175 ERR260139 0.0000000000 -8.1901366    T2D met-

dim(data_in) # 492408      5

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

mean(no_compounds) # 4975.972
sd(no_compounds) #  489.4929

mean(sample_sum_relabun) # 70.46745
sd(sample_sum_relabun) # 2.930408

length(unique(data_in$cpd_id)) # 6839

#-------------------------

# all 3 P < 0.05
#### Forslund T2D-SWE - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (>= 5th percentile)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty72-Hostremoval-EVEN-seqs-5th-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-5th-qty72-v8e.rds")
str(df) # 'data.frame':	492408 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2592 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  8.46e-06 7.74e-06 7.74e-06 7.74e-06 7.74e-06 ...
# $ log10_abun  : num  -5.07 -5.11 -5.11 -5.11 -5.11 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 30
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 42

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 470, p-value = 0.03414
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.0075 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
# D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-5th-qty72-v8e.rds")
str(df) # 'data.frame':	492408 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	360 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd03198" "cpd00224" "cpd00382" "cpd00082" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.0772 0.1288 0.0778 0.2189 0.1046 ...
# $ log10_abun  : num  -1.11 -0.89 -1.11 -0.66 -0.98 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 841, p-value = 0.007768
# alternative hypothesis: true location shift is greater than 0


test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), 0.63 )+
  expand_limits(y = 1.16*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-5th-qty72-v8e.rds")
str(df) # 492408 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	216 obs. of  5 variables:
# $ cpd_id      : chr  "cpd12745" "cpd01554" "cpd01722" "cpd12745" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -8.19 -8.19 -8.19 -8.19 -8.19 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE Rarefied (P5)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	72 obs. of  5 variables:
#   $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.002786 0.000177 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 2 1 1 2 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE Rarefied (P5)" "T2D-SWE Rarefied (P5)" "T2D-SWE Rarefied (P5)" "T2D-SWE Rarefied (P5)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -4.252  -2.555 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.855  -5.855  -5.855  -5.087  -4.252  -2.555 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 487.5, p-value = 0.03503
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-5thperc-even-sequences-T2D-SWE-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------




##########################
########################## T2D-SWE MIN
##########################
##########################


#### T2D Swedish (SWE) cohort - RERUN subset with even sequences

#### Forslund T2D-SWE - w/ Host-removal - keep all samples but rarefy to minimum library size
#-------------------------

#saveRDS(non_host_reads, "non_host_reads.forslund-t2d-swe.rds")
non_host_reads <- readRDS("non_host_reads.forslund-t2d-swe.rds")

hist(non_host_reads);summary(non_host_reads)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 1223102  5572690  7820878  9073574 12868662 22466068 

sum(non_host_reads) # 689591603 = 689,591,603

length(non_host_reads) # 76

min(non_host_reads) # 1223102 = 1,223,102


keep_t2d_swe_list_min <- names(non_host_reads)


writeLines(keep_t2d_swe_list_min, con = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/keep_t2d_swe_list_min.txt")

#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - read in superfocus - fxn potential outputs - RERUN subset with even sequences (minimum library size)
#-------------------------

# SUPER-FOCUS results copied here ...

superfocus_out_dir <- "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min"
list.dirs(superfocus_out_dir)
head( list.dirs(superfocus_out_dir) )

# don't keep 1st directory
( results_dirs <- list.dirs(superfocus_out_dir)[-c(1)] )
length(results_dirs) # 76

head(results_dirs)
# [1] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260139"
# [2] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260140"
# [3] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260144"
# [4] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260147"
# [5] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260151"
# [6] "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260152"

names(results_dirs) <- gsub(pattern = "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_", replacement = "", x = results_dirs)
head(results_dirs)
# ERR260139 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260139" 
# ERR260140 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260140" 
# ERR260144 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260144" 
# ERR260147 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260147" 
# ERR260151 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260151" 
# ERR260152 
# "/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/3_fxn_superfocus_copy_hostremoval_min/superfocus_out_ERR260152" 

sampid <- keep_t2d_swe_list_min

# check identical order
identical(sampid, names(results_dirs)) # FALSE
identical(sort(sampid), sort(names(results_dirs))) # TRUE
length(results_dirs) # 76


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
  
  
  tab$sampid <- this_samp
  names(tab)
  
  #tab <- tab[,c(7,1,2,3,4,6)]
  
  # last column is sampid
  # take average of percentages
  
  #sel.col.percent <- grep(pattern = "_non_host.1.fastq..$", x = names(tab))
  sel.col.percent <- grep(pattern = "_non_host_rarefy_even.1.fastq..$", x = names(tab))
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
dim(sfx.long) # 424919      6
head(sfx.long)
# sampleID                   subsys_L1                    subsys_L2                           subsys_L3                                                      fxn percent_abun
# 2 ERR260268 Amino Acids and Derivatives                            -                 Amino acid racemase                            Alanine_racemase_(EC_5.1.1.1)   0.03601873
# 3 ERR260268 Amino Acids and Derivatives                            -                 Amino acid racemase                   Diaminopimelate_epimerase_(EC_5.1.1.7)   0.01200624
# 4 ERR260268 Amino Acids and Derivatives                            - Creatine and Creatinine Degradation                  Creatinine_amidohydrolase_(EC_3.5.2.10)   0.01200624
# 5 ERR260268 Amino Acids and Derivatives Alanine, serine, and glycine                Alanine biosynthesis                            Alanine_racemase_(EC_5.1.1.1)   0.04802497
# 6 ERR260268 Amino Acids and Derivatives Alanine, serine, and glycine                Alanine biosynthesis Branched-chain_amino_acid_aminotransferase_(EC_2.6.1.42)   0.04722456
# 7 ERR260268 Amino Acids and Derivatives Alanine, serine, and glycine                Alanine biosynthesis                        Cysteine_desulfurase_(EC_2.8.1.7)   0.01460760

sfx.long$full_fxn_tax <- paste0(sfx.long$subsys_L1,"___", sfx.long$subsys_L2,"___", sfx.long$subsys_L3,"___", sfx.long$fxn)

## translate from long to wide format
names(sfx.long)
# "sampleID"     "subsys_L1"    "subsys_L2"    "subsys_L3"    "fxn"          "percent_abun" "full_fxn_tax"

sfx.wide <- dcast(sfx.long, formula = full_fxn_tax ~ sampleID, value.var = "percent_abun")
dim(sfx.wide) #  14820    77

sel.na <- which(is.na(sfx.wide),arr.ind = TRUE)
sfx.wide[sel.na] <- 0

# function taxonomy
full_fxn_names <- sfx.wide$full_fxn_tax

length(full_fxn_names) # 14820
length(unique(full_fxn_names)) # 14820

names(full_fxn_names) <- paste0("fxn_",c(1:length(full_fxn_names)))
head(full_fxn_names)
# fxn_1 
# "Amino Acids and Derivatives___-___Amino acid racemase___4-hydroxyproline_epimerase_(EC_5.1.1.8)" 
# fxn_2 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)" 
# fxn_3 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic" 
# fxn_4 
# "Amino Acids and Derivatives___-___Amino acid racemase___Alanine_racemase_(EC_5.1.1.1)_##_catabolic" 
# fxn_5 
# "Amino Acids and Derivatives___-___Amino acid racemase___Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)" 
# fxn_6 
# "Amino Acids and Derivatives___-___Amino acid racemase___Aspartate_racemase_(EC_5.1.1.13)" 


tax.fxn <- separate(sfx.wide, full_fxn_tax, c("subsys_L1", "subsys_L2", "subsys_L3", "fxn"), sep= "___", remove=TRUE)
# remove sample ids
tax.fxn <- tax.fxn[ ,-which(names(tax.fxn) %in% sampid)]

row.names(tax.fxn) <- names(full_fxn_names)

head(sfx.wide)

names(sfx.wide)
# [1] "full_fxn_tax" "ERR260139"    "ERR260140"    "ERR260144"    "ERR260147"    "ERR260151"    "ERR260152"    "ERR260153"    "ERR260159"    "ERR260161"    "ERR260162"   
# [12] "ERR260163"    "ERR260165"    "ERR260166"    "ERR260167"    "ERR260169"    "ERR260170"    "ERR260171"    "ERR260173"    "ERR260174"    "ERR260175"    "ERR260179"   
# [23] "ERR260180"    "ERR260181"    "ERR260185"    "ERR260186"    "ERR260188"    "ERR260189"    "ERR260190"    "ERR260193"    "ERR260198"    "ERR260199"    "ERR260201"   
# [34] "ERR260203"    "ERR260204"    "ERR260205"    "ERR260206"    "ERR260207"    "ERR260209"    "ERR260210"    "ERR260215"    "ERR260216"    "ERR260217"    "ERR260218"   
# [45] "ERR260221"    "ERR260223"    "ERR260224"    "ERR260225"    "ERR260226"    "ERR260227"    "ERR260230"    "ERR260231"    "ERR260234"    "ERR260241"    "ERR260242"   
# [56] "ERR260243"    "ERR260244"    "ERR260246"    "ERR260250"    "ERR260251"    "ERR260252"    "ERR260253"    "ERR260255"    "ERR260256"    "ERR260258"    "ERR260259"   
# [67] "ERR260260"    "ERR260263"    "ERR260264"    "ERR260265"    "ERR260266"    "ERR260267"    "ERR260268"    "ERR260271"    "ERR260273"    "ERR260276"    "ERR275252"  

#names(sfx.wide) <- gsub(pattern = "-", replacement = "_", x = names(sfx.wide))

identical(as.character(full_fxn_names), sfx.wide$full_fxn_tax) # TRUE

row.names(sfx.wide) <- names(full_fxn_names)
sfx.wide <- sfx.wide[ ,-1]

names(sfx.wide)

head(sampid)
# "ERR260268" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255"

length(sampid) # 76

names(sampid) # NULL - in this case there is NOT an alternative sample name being used

# check alignment of sample IDs and sample names
identical(names(sfx.wide) , as.character(sampid)) # FALSE
identical(sort(names(sfx.wide)) , sort(as.character(sampid))) # TRUE

#NOT RUN THIS TIME
#names(sfx.wide) <- names(sampid)


names(tax.fxn) # "subsys_L1" "subsys_L2" "subsys_L3" "fxn"
dim(tax.fxn) # 14820     4

length(unique(tax.fxn$subsys_L1)) # 35
length(unique(tax.fxn$subsys_L2)) # 176
length(unique(tax.fxn$subsys_L3)) # 990
length(unique(tax.fxn$fxn)) # 8044


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - functions - get into Phyloseq object - RERUN subset with even sequences (minimum library size)
#-------------------------

# sfx.wide - is equiv to OTU table

# tax.fxn - is equiv to TAX table

# meta - is equiv to sample table

## Create 'taxonomyTable'
#  tax_table - Works on any character matrix. 
#  The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
tax.m <- as.matrix( tax.fxn )
dim(tax.m) # 14820     4

TAX <- tax_table( tax.m )


## Create 'otuTable'
#  otu_table - Works on any numeric matrix. 
#  You must also specify if the species are rows or columns
otu.m <- as.matrix( sfx.wide )
dim(otu.m)
# 14820    76

OTU <- otu_table(otu.m, taxa_are_rows = TRUE)


## Create a phyloseq object, merging OTU & TAX tables
phy = phyloseq(OTU, TAX)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14820 taxa and 76 samples ]
# tax_table()   Taxonomy Table:    [ 14820 taxa by 4 taxonomic ranks ]

sample_names(phy)
# [1] "ERR260139" "ERR260140" "ERR260144" "ERR260147" "ERR260151" "ERR260152" "ERR260153" "ERR260159" "ERR260161" "ERR260162" "ERR260163" "ERR260165" "ERR260166" "ERR260167"
# [15] "ERR260169" "ERR260170" "ERR260171" "ERR260173" "ERR260174" "ERR260175" "ERR260179" "ERR260180" "ERR260181" "ERR260185" "ERR260186" "ERR260188" "ERR260189" "ERR260190"
# [29] "ERR260193" "ERR260198" "ERR260199" "ERR260201" "ERR260203" "ERR260204" "ERR260205" "ERR260206" "ERR260207" "ERR260209" "ERR260210" "ERR260215" "ERR260216" "ERR260217"
# [43] "ERR260218" "ERR260221" "ERR260223" "ERR260224" "ERR260225" "ERR260226" "ERR260227" "ERR260230" "ERR260231" "ERR260234" "ERR260241" "ERR260242" "ERR260243" "ERR260244"
# [57] "ERR260246" "ERR260250" "ERR260251" "ERR260252" "ERR260253" "ERR260255" "ERR260256" "ERR260258" "ERR260259" "ERR260260" "ERR260263" "ERR260264" "ERR260265" "ERR260266"
# [71] "ERR260267" "ERR260268" "ERR260271" "ERR260273" "ERR260276" "ERR275252"


### Now Add sample data to phyloseq object
# sample_data - Works on any data.frame. The rownames must match the sample names in
# the otu_table if you plan to combine them as a phyloseq-object

# reuse subset of previous fxn phyloseq object
temp <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-v8d.RDS")
temp <- prune_samples(samples = sample_names(phy), x = temp)

df.samp <- as(temp@sam_data, "data.frame")

head(df.samp)
# Sample Country.subset         Status      Bases       Run group_label      age non_host_reads fxn_sum_counts
# ERR260139 NG-5636_334            SWE T2D metformin- 2036676514 ERR260139    T2D met- 70.25205        5248535         198299
# ERR260140 NG-5636_344            SWE T2D metformin- 1935856900 ERR260140    T2D met- 70.15342        5378909         255775
# ERR260144 NG-5636_353            SWE T2D metformin- 2483902494 ERR260144    T2D met- 69.57534        7906030         268694
# ERR260147 NG-5636_365            SWE        ND CTRL 2821768300 ERR260147      Normal 71.39452        6729275         374375
# ERR260151 NG-5636_378            SWE T2D metformin- 2630431274 ERR260151    T2D met- 71.56712        7922852         200815
# ERR260152 NG-5636_380            SWE T2D metformin- 1813559434 ERR260152    T2D met- 71.24384        4516029         154717

# remove columns: 'Bases', non_host_reads, fxn_sum_counts, as not applicable to this version based on rarefied sequences
dim(df.samp) #  76  9
sel <- which(names(df.samp) %in% c("Bases","non_host_reads","fxn_sum_counts"))
df.samp <- df.samp[ ,-sel]
head(df.samp)

# reorder to align with phy object
df.samp2 <- df.samp[ sample_names(phy), ]
identical(row.names(df.samp2), sample_names(phy)) # TRUE

SAMP <- sample_data(df.samp2)


### Combine SAMPDATA into phyloseq object
phy <- merge_phyloseq(phy, SAMP)
phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14820 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 14820 taxa by 4 taxonomic ranks ]

head(taxa_names(phy))
# "fxn_1" "fxn_2" "fxn_3" "fxn_4" "fxn_5" "fxn_6"

head(phy@tax_table)
# Taxonomy Table:     [6 taxa by 4 taxonomic ranks]:
#   subsys_L1                     subsys_L2 subsys_L3             fxn                                                                                               
# fxn_1 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "4-hydroxyproline_epimerase_(EC_5.1.1.8)"                                                         
# fxn_2 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)"                                                                   
# fxn_3 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_biosynthetic"                                                   
# fxn_4 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Alanine_racemase_(EC_5.1.1.1)_##_catabolic"                                                      
# fxn_5 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Arginine_racemase_(EC_5.1.1.9)_@_Lysine_racemase_(EC_5.1.1.5)_@_Ornithine_racemase_(EC_5.1.1.12)"
# fxn_6 "Amino Acids and Derivatives" "-"       "Amino acid racemase" "Aspartate_racemase_(EC_5.1.1.13)"                         

getwd()  # "/Users/lidd0026/WORKSPACE/PROJ/PCaN-NZ/nz-city-resto/modelling/R"

table(phy@sam_data$group_label)
# T2D met-   Normal 
# 33       43 

saveRDS(object = phy, file = "phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS")

#phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS")

# get stats?

head(phy@otu_table)
fxns <- as.data.frame( phy@otu_table )
NonZeroFxns <- apply( fxns , 2,function(x) length(which(x > 0)) )
length(NonZeroFxns) # 76
NonZeroFxns

mean(NonZeroFxns) # 5591.039
sd(NonZeroFxns) # 1169.659


#-------------------------

#### Forslund T2D-SWE - w/ Host removal - COPY of R code to run CPP steps on HPC - RERUN subset with even sequences (minimum library size)
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
# # For study - Forslund et al T2D-SWE rarefied sequences - minimum library size
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
# message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min'")
# workdir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min"
# message("\nsetwd(workdir)")
# setwd(workdir)
# message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min/working'")
# temp_dir <- "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min/working"
# 
# message("\nthis_study <- '-t2d-swe-rarefied-min-pawsey'")
# this_study <- "-t2d-swe-rarefied-min-pawsey"
# message("\nphy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS')")
# phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS")
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

#### Forslund T2D-SWE - w/ Host-removal - COPY of OUTOUTS from R code after running CPP steps on HPC - RERUN subset with even sequences (minimum library size)
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
# workdir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min'
# 
# setwd(workdir)
# 
# temp_dir <- '/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min/working'
# 
# this_study <- '-t2d-swe-rarefied-min-pawsey'
# 
# phy <- readRDS('phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS')
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
# [1] 14820     4
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# ...
# 
# 
# 
# [[14818]]
# NULL
# 
# [[14819]]
# NULL
# 
# [[14820]]
# NULL
# 
# 
# ## assemble results
# 
# (num_results_files <- dim(df.tax)[1])
# [1] 14820
# [1] "added df 1 of 14820"
# [1] "added df 2 of 14820"
# [1] "added df 3 of 14820"
# ...
# 
# 
# [1] "added df 14818 of 14820"
# [1] "added df 14819 of 14820"
# [1] "added df 14820 of 14820"
# 
# str(df.out)
# 'data.frame':	388231 obs. of  8 variables:
#   $ superfocus_fxn: chr  NA "fxn_1" "fxn_1" "fxn_2" ...
# $ f             : int  NA 1 1 1 1 1 1 1 1 1 ...
# $ f__in         : chr  NA "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id        : chr  NA "rxn02360" "rxn02360" "rxn00283" ...
# $ cpd_id        : chr  NA "cpd00851" "cpd02175" "cpd00035" ...
# $ cpd_name      : chr  NA "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" ...
# $ cpd_form      : chr  NA "C5H9NO3" "C5H9NO3" "C3H7NO2" ...
# $ cpd_molar_prop: num  NA 1 1 1 1 1 1 1 1 1 ...
# 
# head(df.out)
# superfocus_fxn  f                                   f__in   rxn_id   cpd_id
# 1           <NA> NA                                    <NA>     <NA>     <NA>
#   2          fxn_1  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_1  1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_2  1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# cpd_name cpd_form cpd_molar_prop
# 1                      <NA>     <NA>             NA
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1
# 4                 L-Alanine  C3H7NO2              1
# 5                 D-Alanine  C3H7NO2              1
# 6                 L-Alanine  C3H7NO2              1
# 
# dim(df.out)
# [1] 388230      8
# 
# ## normalise molar_prop to cpd_relabun so total of 1 per superfocus function
# 
# length(unique(df.out$superfocus_fxn))
# [1] 8318
# 
# phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 14820 taxa and 76 samples ]
# sample_data() Sample Data:       [ 76 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 14820 taxa by 4 taxonomic ranks ]
# 
# % of functions represented - with compound information
# [1] 56.12686
# [1] "completed 1"
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 8316"
# [1] "completed 8317"
# [1] "completed 8318"
# 
# sum(df.out$cpd_molar_prop_norm)
# [1] 8318
# 
# sample_sums(phy)
# ERR260139 ERR260140 ERR260144 ERR260147 ERR260151 ERR260152 ERR260153 ERR260159 
# 100       100       100       100       100       100       100       100 
# ERR260161 ERR260162 ERR260163 ERR260165 ERR260166 ERR260167 ERR260169 ERR260170 
# 100       100       100       100       100       100       100       100 
# ERR260171 ERR260173 ERR260174 ERR260175 ERR260179 ERR260180 ERR260181 ERR260185 
# 100       100       100       100       100       100       100       100 
# ERR260186 ERR260188 ERR260189 ERR260190 ERR260193 ERR260198 ERR260199 ERR260201 
# 100       100       100       100       100       100       100       100 
# ERR260203 ERR260204 ERR260205 ERR260206 ERR260207 ERR260209 ERR260210 ERR260215 
# 100       100       100       100       100       100       100       100 
# ERR260216 ERR260217 ERR260218 ERR260221 ERR260223 ERR260224 ERR260225 ERR260226 
# 100       100       100       100       100       100       100       100 
# ERR260227 ERR260230 ERR260231 ERR260234 ERR260241 ERR260242 ERR260243 ERR260244 
# 100       100       100       100       100       100       100       100 
# ERR260246 ERR260250 ERR260251 ERR260252 ERR260253 ERR260255 ERR260256 ERR260258 
# 100       100       100       100       100       100       100       100 
# ERR260259 ERR260260 ERR260263 ERR260264 ERR260265 ERR260266 ERR260267 ERR260268 
# 100       100       100       100       100       100       100       100 
# ERR260271 ERR260273 ERR260276 ERR275252 
# 100       100       100       100 
# 
# getwd()
# [1] "/scratch/pawsey1216/cliddicoat/ft2d_swe/cpp_analysis_min"
# 
# ### 2) get cpd rel abun per sample
# 
# # # # # # # # # # #
# 
# dim(df.OTU)
# [1] 14820    76
# [[1]]
# NULL
# 
# [[2]]
# NULL
# 
# [[3]]
# NULL
# 
# ...
# 
# 
# 
# [[75]]
# NULL
# 
# [[76]]
# NULL
# 
# 
# ## assemble results
# superfocus_fxn f                                   f__in   rxn_id   cpd_id
# 2          fxn_1 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd00851
# 3          fxn_1 1 4-hydroxyproline epimerase (EC 5.1.1.8) rxn02360 cpd02175
# 4          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00035
# 5          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn00283 cpd00117
# 6          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00035
# 7          fxn_2 1           Alanine racemase (EC 5.1.1.1) rxn19085 cpd00117
# cpd_name cpd_form cpd_molar_prop cpd_molar_prop_norm
# 2 trans-4-Hydroxy-L-proline  C5H9NO3              1           0.5000000
# 3   cis-4-Hydroxy-D-proline  C5H9NO3              1           0.5000000
# 4                 L-Alanine  C3H7NO2              1           0.1666667
# 5                 D-Alanine  C3H7NO2              1           0.1666667
# 6                 L-Alanine  C3H7NO2              1           0.1666667
# 7                 D-Alanine  C3H7NO2              1           0.1666667
# sample cpd_rel_abun_norm
# 2 ERR260139       0.000000000
# 3 ERR260139       0.000000000
# 4 ERR260139       0.001681027
# 5 ERR260139       0.001681027
# 6 ERR260139       0.001681027
# 7 ERR260139       0.001681027
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 74"
# [1] "completed 75"
# [1] "completed 76"
# 
# str(dat)
# 'data.frame':	29505480 obs. of  11 variables:
#   $ superfocus_fxn     : chr  "fxn_1" "fxn_1" "fxn_2" "fxn_2" ...
# $ f                  : int  1 1 1 1 1 1 1 1 1 1 ...
# $ f__in              : chr  "4-hydroxyproline epimerase (EC 5.1.1.8)" "4-hydroxyproline epimerase (EC 5.1.1.8)" "Alanine racemase (EC 5.1.1.1)" "Alanine racemase (EC 5.1.1.1)" ...
# $ rxn_id             : chr  "rxn02360" "rxn02360" "rxn00283" "rxn00283" ...
# $ cpd_id             : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ cpd_name           : chr  "trans-4-Hydroxy-L-proline" "cis-4-Hydroxy-D-proline" "L-Alanine" "D-Alanine" ...
# $ cpd_form           : chr  "C5H9NO3" "C5H9NO3" "C3H7NO2" "C3H7NO2" ...
# $ cpd_molar_prop     : num  1 1 1 1 1 1 1 1 1 1 ...
# $ cpd_molar_prop_norm: num  0.5 0.5 0.167 0.167 0.167 ...
# $ sample             : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun_norm  : num  0 0 0.00168 0.00168 0.00168 ...
# 
# sum(dat$cpd_rel_abun_norm)
# [1] 5357.555
# 
# average functional relative abundance per sample
# 
# sum(dat$cpd_rel_abun_norm)/nsamples(phy)
# [1] 70.49415
# 
# names(dat)
# [1] "superfocus_fxn"      "f"                   "f__in"              
# [4] "rxn_id"              "cpd_id"              "cpd_name"           
# [7] "cpd_form"            "cpd_molar_prop"      "cpd_molar_prop_norm"
# [10] "sample"              "cpd_rel_abun_norm"  
# 
# length(unique(dat$cpd_id))
# [1] 6785
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
# ...
# 
# 
# 
# [[75]]
# NULL
# 
# [[76]]
# NULL
# 
# 
# ## assemble results
# cpd_id    sample cpd_rel_abun
# 1 cpd00851 ERR260139  0.001062626
# 2 cpd02175 ERR260139  0.000000000
# 3 cpd00035 ERR260139  0.261933310
# 4 cpd00117 ERR260139  0.225371027
# 5 cpd00051 ERR260139  0.108537085
# 6 cpd00586 ERR260139  0.000000000
# [1] "completed 2"
# [1] "completed 3"
# ...
# 
# 
# [1] "completed 74"
# [1] "completed 75"
# [1] "completed 76"
# 
# str(dat.cpd.collate)
# 'data.frame':	515660 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.00106 0 0.26193 0.22537 0.10854 ...
# 
# sum(dat.cpd.collate$cpd_rel_abun)
# [1] 5357.555
# 
# sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))
# [1] 70.49415
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries
# [CRAYBLAS_WARNING] Application linked against multiple cray-libsci libraries


#-------------------------

#### Forslund T2D-SWE - w/ Host-removal - continue CPP analysis - RERUN subset with even sequences (minimum library size)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS")

# copy output file from HPC
dat.cpd.collate <- readRDS("/Users/lidd0026/WORKSPACE/PROJ/Gut-and-soil/modelling/PawseyHPCSupp/forslund-t2d-swe-EVEN-sequences/cpp_analysis_min/dat.cpd.collate-all-samps-cpp3d--t2d-swe-rarefied-min-pawsey.rds")

hist(dat.cpd.collate$cpd_rel_abun); summary(dat.cpd.collate$cpd_rel_abun)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.000000 0.000172 0.010390 0.001668 7.167470 

length(unique(dat.cpd.collate$cpd_id)) # 6785
length(unique(dat.cpd.collate$sample)) # 76
str(dat.cpd.collate)
# 'data.frame':	515660 obs. of  3 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.00106 0 0.26193 0.22537 0.10854 ...
6785*76 # 515660

hist(log10(dat.cpd.collate$cpd_rel_abun)); summary(log10(dat.cpd.collate$cpd_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf -3.7646    -Inf -2.7777  0.8554 

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
# -7.6831 -7.6831 -3.7646 -4.5425 -2.7777  0.8554 

# make group variable from sample name

dat.cpd.collate$group_label <- NA

# from above
df.samp <- as(phy@sam_data, "data.frame")

identical( phy@sam_data$Run , df.samp$Run ) # TRUE
identical( sample_names(phy), df.samp$Run ) # TRUE
unique(df.samp$group_label)
# [1] T2D met- Normal  
# Levels: T2D met- < Normal

#for (i in 1:length(sample_names(phy))) {
for (i in 1:length( df.samp$Run )) {
  #i<-1
  #this_samp <- sample_names(phy)[i]
  this_samp <- df.samp$Run[i]
  sel <- which(dat.cpd.collate$sample == this_samp)
  #dat.cpd.collate$group[sel] <- phy@sam_data$age[i]
  dat.cpd.collate$group_label[sel] <- as.character( df.samp$group_label[i] )
  print(paste0("completed ", i))
}

unique(dat.cpd.collate$group_label) # "T2D met-" "Normal"
dat.cpd.collate$group_label <- factor(dat.cpd.collate$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)
head(dat.cpd.collate)

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-min-qty76-v8e.rds" )
#dat.cpd.collate <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-min-qty76-v8e.rds")

str(dat.cpd.collate)
# 'data.frame':	515660 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd00851" "cpd02175" "cpd00035" "cpd00117" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.00106 0 0.26193 0.22537 0.10854 ...
# $ log10_abun  : num  -2.974 -7.683 -0.582 -0.647 -0.964 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...


## CPP stats ?

data_in <- dat.cpd.collate

head(data_in)
# cpd_id    sample cpd_rel_abun log10_abun group_label
# 1 cpd00851 ERR260139  0.001062626 -2.9736194    T2D met-
#   2 cpd02175 ERR260139  0.000000000 -7.6831064    T2D met-
#   3 cpd00035 ERR260139  0.261933310 -0.5818093    T2D met-
#   4 cpd00117 ERR260139  0.225371027 -0.6471019    T2D met-
#   5 cpd00051 ERR260139  0.108537085 -0.9644218    T2D met-
#   6 cpd00586 ERR260139  0.000000000 -7.6831064    T2D met-

dim(data_in) # 515660      5

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

mean(no_compounds) # 4706.724
sd(no_compounds) #  491.5195

mean(sample_sum_relabun) # 70.49415
sd(sample_sum_relabun) # 2.960926

length(unique(data_in$cpd_id)) # 6785

#-------------------------

# all 3 have P < 0.05
#### Forslund T2D-SWE - w/ Host-removal - check for robustness of key signals using RERUN subset with even sequences (minimum library size)
#-------------------------

phy <- readRDS("phy-phyloseq-fxn-Forslund-SWE-T2D-qty76-Hostremoval-EVEN-seqs-min-v8e.RDS")
df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-min-qty76-v8e.rds")
str(df) # 'data.frame':	515660 obs. of  5 variables:


## T2D-SWE - BCFA-ACPs

sel <- which(df$cpd_id %in% new_bcfa)
df <- df[sel, ]
length(unique(df$cpd_id)) # 36

str(df)
# 'data.frame':	2736 obs. of  5 variables:
# $ cpd_id      : chr  "cpd11472" "cpd11475" "cpd11465" "cpd11469" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0 0 0 0 0 0 0 0 0 0 ...
# $ log10_abun  : num  -7.68 -7.68 -7.68 -7.68 -7.68 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "BCFA-ACPs"
res$dataset <- "T2D-SWE Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ] # 33
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ] # 43

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 533, p-value = 0.03256
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), max(res$sum_rel_abun) + 0.0008 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(d)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-SWE-BCFA-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")




## T2D-SWE - Sugars
#  D-Fructose = cpd00082 ; L-Arabinose = cpd00224 ; Melibiose = cpd03198 ; 6-Phosphosucrose = cpd01693 ; Melitose (Raffinose) = cpd00382

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-min-qty76-v8e.rds")
str(df) # 'data.frame':	515660 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd00082", "cpd00224", "cpd03198", "cpd01693", "cpd00382"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 5

str(df)
# 'data.frame':	380 obs. of  5 variables:
#   $ cpd_id      : chr  "cpd03198" "cpd00224" "cpd00382" "cpd00082" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260139" ...
# $ cpd_rel_abun: num  0.0786 0.1302 0.0791 0.2292 0.1096 ...
# $ log10_abun  : num  -1.104 -0.885 -1.102 -0.64 -0.96 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 1 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Sugars"
res$dataset <- "T2D-SWE Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)

x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "greater" ,  paired = FALSE) # 
wmw.test
# Wilcoxon rank sum exact test
# data:  x and y
# W = 963, p-value = 0.003744
# alternative hypothesis: true location shift is greater than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,4))

p <- ggplot(data = res, aes(x = group_label, y = sum_rel_abun) )+
  #ylim( min(res$sum_rel_abun), max(res$sum_rel_abun) + 0.07 )+
  expand_limits(y = 1.1*max(res$sum_rel_abun))+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("Summed CPP (%)")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "right", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(e)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-SWE-Sugars-v8h.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")


## T2D-SWE - Lignin\n& precursors
# Lignin = cpd12745 ; Sinapyl alcohol = cpd01554 ; p-Coumaryl alcohol = cpd01722

df <- readRDS("dat.cpd.collate-all-samps-cpp3d--forslund-t2d-swe-hostremoval-ExtraData-EVEN-seqs-min-qty76-v8e.rds")
str(df) # 515660 obs. of  5 variables:

sel <- which(df$cpd_id %in% c( "cpd12745", "cpd01554", "cpd01722"))
df <- df[sel, ]
length(unique(df$cpd_id)) # 3

str(df)
# 'data.frame':	228 obs. of  5 variables:
# $ cpd_id      : chr  "cpd01554" "cpd01722" "cpd12745" "cpd01554" ...
# $ sample      : chr  "ERR260139" "ERR260139" "ERR260139" "ERR260140" ...
# $ cpd_rel_abun: num  0 0 0 0 0 ...
# $ log10_abun  : num  -7.68 -7.68 -7.68 -7.68 -7.68 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 1 1 1 1 1 1 2 ...

#df$group_label <- df$group

res <- data.frame(sample = unique(df$sample), sum_rel_abun = NA, group_label = NA )

for (i in 1:length(unique(df$sample))) {
  #i<-1
  this_samp <- res$sample[i]
  subsel <- which(df$sample == this_samp)
  res$sum_rel_abun[i] <- sum(df$cpd_rel_abun[subsel])
  res$group_label[i] <- as.character(unique(df$group_label[subsel]))
  
  print(paste0("completed ",i))
}

res$cpd_group <- "Lignin & precursors"
res$dataset <- "T2D-SWE Rarefied (Min)"

unique(res$group_label) # "T2D met-" "Normal"  
res$group_label <- factor(res$group_label, levels = c("T2D met-", "Normal"), ordered = TRUE)

str(res)
# 'data.frame':	76 obs. of  5 variables:
#   $ sample      : chr  "ERR260139" "ERR260140" "ERR260144" "ERR260147" ...
# $ sum_rel_abun: num  0 0 0 0.00249 0 ...
# $ group_label : Ord.factor w/ 2 levels "T2D met-"<"Normal": 1 1 1 2 1 1 2 1 1 1 ...
# $ cpd_group   : chr  "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" "Lignin & precursors" ...
# $ dataset     : chr  "T2D-SWE Rarefied (Min)" "T2D-SWE Rarefied (Min)" "T2D-SWE Rarefied (Min)" "T2D-SWE Rarefied (Min)" ...

# use log10 of summed rel abun

hist(log10(res$sum_rel_abun)); summary(log10(res$sum_rel_abun))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -Inf    -Inf    -Inf    -Inf  -4.120  -2.604 

# log10 abun
res$log10_sum_rel_abun <- res$sum_rel_abun
# set zero-replacement value at 1/2 smallest non-zero value of that group
subsel.zero <- which(res$log10_sum_rel_abun == 0) #
if (length(subsel.zero) > 0) {
  zero_replace <- 0.5*min(res$log10_sum_rel_abun[ -subsel.zero ])
  res$log10_sum_rel_abun[ subsel.zero ] <- zero_replace
}
res$log10_sum_rel_abun <- log10(res$log10_sum_rel_abun)

hist(res$log10_sum_rel_abun); summary( res$log10_sum_rel_abun )
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -5.658  -5.658  -5.658  -5.018  -4.120  -2.604 

#x <- res$sum_rel_abun[ which(res$group_label == "T2D met-") ]
#y <- res$sum_rel_abun[ which(res$group_label == "Normal") ]
x <- res$log10_sum_rel_abun[ which(res$group_label == "T2D met-") ]
y <- res$log10_sum_rel_abun[ which(res$group_label == "Normal") ]

wmw.test <- wilcox.test(x, y, alternative = "less" ,  paired = FALSE) # Results are same for Summed CPP% and log10(Summed CPP%)
wmw.test
# Wilcoxon rank sum test with continuity correction
# data:  x and y
# W = 567.5, p-value = 0.03977
# alternative hypothesis: true location shift is less than 0

test_result <- paste0(unique(res$dataset),": ",unique(res$cpd_group),"\n",
                      #"T2D Met- vs Normal (SWE) Rarefied\n",
                      "Wilcoxon-Mann-Whitney\nW = ",round(wmw.test$statistic,0),", P = ",round(wmw.test$p.value,3))

p <- ggplot(data = res, aes(x = group_label, y = log10_sum_rel_abun) )+ # y = sum_rel_abun
  ylim( min(res$log10_sum_rel_abun), -2.3 )+
  geom_violin()+
  geom_boxplot(width = 0.2, alpha = 0.3)+
  geom_jitter(width = 0.1, height = 0, alpha = 0.3)+
  xlab("Diagnosis")+ ylab("log10(Summed CPP (%))")+
  theme_bw()+
  annotate(geom="text_npc", npcx = "left", npcy = "top", label = test_result, size = 2.75 , lineheight = 0.85)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = rel(1.1)),
    #axis.text.x = element_text(size = rel(0.9), angle = 15, hjust=1, vjust=1),
    #plot.title = element_text(hjust = 0.5, size = rel(1)),
    axis.title = element_text(size = rel(0.9))
  )

p

grid.text(label = "(f)", x = unit(0.04, "npc") , y = unit(0.96,"npc"), gp=gpar(fontsize=13, fontface="bold") )
dev.print(tiff, file = paste0(workdir,"/plots/","Rarefied-min-even-sequences-T2D-SWE-Lignin&precursors-v8e.tiff"), width = 8, height = 8, units = "cm", res=600, compression="lzw",type="cairo")



#-------------------------


