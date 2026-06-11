############
#
# R script for cpp3d
# - build reaction search in parallel - get_reactions & compounds
# - get cpd rel abun per sample
# - collate_compounds
#
# For study - AMI disturbed vs natural soils
# Craig Liddicoat - Flinders University
# Running on Pawsey Setonix
############

# Add a new path
.libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
            "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))

R.Version()

# load packages
#library(readxl); packageVersion("readxl")
library(parallel); packageVersion("parallel")
library(doParallel); packageVersion("doParallel")
library(dplyr); packageVersion("dplyr")
library(stringr); packageVersion("stringr")
library(phyloseq); packageVersion("phyloseq") # '1.44.0'

message("\n# establish folders and input files")

message("\nworkdir <- '/scratch/pawsey1216/cliddicoat/ami_2025/cpp_analysis'")
workdir <- "/scratch/pawsey1216/cliddicoat/ami_2025/cpp_analysis"
message("\nsetwd(workdir)")
setwd(workdir)
message("\ntemp_dir <- '/scratch/pawsey1216/cliddicoat/ami_2025/cpp_analysis/working'")
temp_dir <- "/scratch/pawsey1216/cliddicoat/ami_2025/cpp_analysis/working"

message("\nthis_study <- '-ami-dist-vs-nat-pawsey'")
this_study <- "-ami-dist-vs-nat-pawsey"
message("\nphy <- readRDS('phy-phyloseq-fxn-ami-dist-vs-nat-v8c.RDS')")
phy <- readRDS("phy-phyloseq-fxn-ami-dist-vs-nat-v8c.RDS")


subsys.lut <- readRDS("subsys.lut.RDS")
rxns.lut <- readRDS("rxns.lut.RDS")
rxn_pathways.lut <- readRDS("rxn_pathways.lut.RDS")
compounds.lut <- readRDS("compounds.lut.RDS")



message("\n### 1) build reaction search in parallel - get_reactions & compounds")
message("\n# # # # # # # # # #")
message("\ndf.tax <- as.data.frame(phy@tax_table)")
df.tax <- as.data.frame(phy@tax_table)
message("\nhead(row.names(df.tax))")
head(row.names(df.tax))
message("\ndim(df.tax)")
dim(df.tax)


get_rxns_and_compounds_indiv <- function( df.tax, subsys.lut, rxns.lut, rxn_pathways.lut ) {
  
  rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
  rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
 
  sub1 <- df.tax$subsys_L1[i]
  sub2 <- df.tax$subsys_L2[i]
  sub3 <- df.tax$subsys_L3[i]
  
  fxn.temp <- df.tax$fxn[i]
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  
  # store results corresponding to each Superfocus row
  fxn.list <- list()
  fxn.list[[ fxn.superfocus.rowlabel  ]] <- list()
  
  # check for multiple functions/reactions?
  flag1 <- grepl(pattern = "_/_|/", x = fxn.temp)
  flag2 <- grepl(pattern = "_@_", x = fxn.temp)
  if (!any(flag1,flag2)==TRUE) {
    # no multiples
    fxns <- fxn.temp
  } else if (flag1==TRUE) {
    fxns <- unlist( strsplit(fxn.temp, split = "_/_") )  ###### WHAT ABOUT SPLIT FOR "/" WITHOUT UNDERSCORES ??
  } else {
    fxns <- unlist( strsplit(fxn.temp, split = "_@_") )
  }
  # remove underscores
  ( fxns <- gsub(pattern = "_", replacement = " ", x = fxns) )
  
  # process each fxn & store attributes
  df.fxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=1:length(fxns),`f__in`=fxns, matching_method=NA, rxns=NA)
  
  # Identify '/' separators with no '_'  ??
  
  for (f in 1:length(fxns)) {  # this accounts for multiple functions/reactions reported in Superfocus outputs
    #f<-1
    #f<-2
    f.in <- fxns[f]
    
    # these concatenated expressions will be used to look for exact match using hierarchy in ModelSEED Subsystem table
    full_hier_target <- paste0(sub1,"__",sub2,"__",sub3,"__",f.in)
    full_hier_list <- paste0(subsys.lut$Class,"__",subsys.lut$Subclass,"__",gsub("_"," ",subsys.lut$Name),"__",subsys.lut$Role)
    
    ## data cleaning
    
    # trim off '_#' and '_##' tags
    trim_nchar <- str_locate(string = f.in, pattern = " # | ## ")[1]
    if (!is.na(trim_nchar) & length(trim_nchar)==1) {
      f.in <- substring(text = f.in , first = 1, last = trim_nchar-1)
    }
    
    # Eliminate unwanted parsing of regular expressions: '[', ']','***', '(', ')'
    f.in <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\} ", replacement ="." , x = f.in) # used later
    
    #rxns.lut$name <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$name) # used later
    #rxns.lut$aliases <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = rxns.lut$aliases) # used later
    
    full_hier_target <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_target)
    full_hier_list <- gsub(pattern = "\\[|\\]|\\*+|\\(|\\)|\\{|\\}", replacement ="." , x = full_hier_list)
    
    sel.rx <- grep(pattern = full_hier_target, x = full_hier_list)
    
    ## ALTERNATIVE #1 == FULL HIERACHICAL MATCH
    if (length(sel.rx)>=1) {
      df.fxns$matching_method[f] <- "Exact hierachy match"
      df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
      
    } else if (str_detect(string = fxns[f], pattern = " \\(EC ")) {  ## ALTERNATIVE #2 == MATCHING ECs
      # search by EC id if present
      
      f.in <- fxns[f] # this goes back to string with brackets for EC
      ## LOOK FOR MULTIPLE ECs ?
      
      how_many_ECs <- str_count(string = f.in, pattern = "\\(EC.*?\\)")
      
      ECs <- as.character( str_extract_all(string = f.in, pattern = "\\(EC.*?\\)", simplify = TRUE) )
      #class(ECs)
      ECs <- gsub(pattern = "\\(EC |\\)", replacement = "", x = ECs)
      ECs.collapse <- paste0(ECs, collapse = "|")
      
      sel.rx <- which(rxns.lut$ec_numbers == ECs.collapse)
      
      if (length(how_many_ECs)==0 | length(ECs)==0) {
        # there was a glitch, database typo, or some error in identifying the EC number
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      } else if (length(sel.rx)>=1) {
        # combined EC hits identified
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(which(rxns.lut$ec_numbers %in% ECs)) >=1) {
        # treat EC hits individually
        sel.rx <- which(rxns.lut$ec_numbers %in% ECs) # look 1st where ECs are exact matches for EC numbers in Reactions lookup table
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if (length(grep(pattern = ECs, x = rxns.lut$ec_numbers)) >=1) {
        # this allows EC to be part of a combination of EC numbers that are listed in Reactions lookup table
        sel.rx <- grep(pattern = ECs, x = rxns.lut$ec_numbers)
        
        df.fxns$matching_method[f] <- "EC number"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else {
        # it had an EC number but couldn't find a match in the EC numbers listed in Reaction lookup table
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      # END EC matching
      
      
    } else {  ## ALTERNATIVE 3 == FXN NAME MATCHING
      ## otherwise attempt to match function name - a) first look for exact matches   ########## then b) closest match above a threshold
      # 1. 'reactions' table by name: rxns.lut$name
      # 2. 'reactions' table by aliases: rxns.lut$aliases
      # 3. 'Model SEED Subsystems' table by Role: subsys.lut$Role
      # 4. 'Unique_ModelSEED_Reaction_Pathways' table by External ID: rxn_pathways.lut$External_rxn_name
      
      if ( length( grep(pattern = f.in, x = rxns.lut$name) )>=1 ) {
        # 1a - exact match - rxns.lut$name
        sel.rx <- grep(pattern = f.in, x = rxns.lut$name)
        #rxns.lut$name[sel.rx]
        df.fxns$matching_method[f] <- "Matched Reactions name"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxns.lut$aliases) )>=1 ) {
        # 2a - exact match - rxns.lut$aliases
        sel.rx <- grep(pattern = f.in, x = rxns.lut$aliases)
        #rxns.lut$aliases[sel.rx]
        #rxns.lut$name[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Reactions aliases"
        df.fxns$rxns[f] <- paste0( unique(rxns.lut$id[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = subsys.lut$Role) )>=1 ) {
        # 3a - exact match - subsys.lut$Role
        sel.rx <- grep(pattern = f.in, x = subsys.lut$Role)
        #subsys.lut$Role[sel.rx]
        #subsys.lut$Reaction[sel.rx]
        
        df.fxns$matching_method[f] <- "Matched Subsytem role"
        df.fxns$rxns[f] <- paste0( unique(subsys.lut$Reaction[sel.rx]), collapse = ";")
        
      } else if ( length( grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name) )>=1 ) {
        # 4a - exact match - rxn_pathways.lut$External_rxn_name
        sel.rx <- grep(pattern = f.in, x = rxn_pathways.lut$External_rxn_name)
        
        df.fxns$matching_method[f] <- "Matched ModelSEED Reaction pathways"
        df.fxns$rxns[f] <- paste0( unique(rxn_pathways.lut$rxn_id[sel.rx]), collapse = ";")
        
        
      } else {
        df.fxns$matching_method[f] <- "No match found"
        df.fxns$rxns[f] <- NA
        
      }
      
      ## DON'T RUN PARTIAL MATCHING AT THIS STAGE
      
      
    } # END function - reaction search
    
    #fxn.list[[ fxn.superfocus.rowlabel  ]][[ f ]][[ "fxns" ]] <- df.fxns
    
    #print(paste0("completed fxn ", f))
    
    
    ## now investigate these reactions ...
    # Reactions lookup table: 
    # - "equation": Definition of reaction expressed using compound IDs and after protonation
    # Compounds lookup table:
    # - "formula": Standard chemical format (using Hill system) in protonated form to match reported charge
    #df.fxns
    
    
    #if (df.fxns$matching_method == "No match found") {
    if (df.fxns$rxns[f] == "" | is.na(df.fxns$rxns[f])) {
      
      df.Rxns <- NA
      df.Compounds <- NA
      
    } else { # reaction(s) were identified
      
      # consider reactions for this f.in only (possibly > 1 f.in per Superfocus row)
      f.in.rxns <- unique(unlist(str_split(string = df.fxns$rxns[f], pattern = ";")))
      
      df.Rxns <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel, f=f, f__in=f.in,rxn_id= f.in.rxns,
                            rxn_name=NA, rxn_eqn=NA, rxn_defn=NA,compds=NA,compd_coef=NA, chem_formx=NA )
                  
      for (r in 1:dim(df.Rxns)[1]) {
        #r<-1
        #this_rxn <- "rxn00004"
        this_rxn <- df.Rxns$rxn_id[r]
        sel <- which(rxns.lut$id == this_rxn)
        ( df.Rxns$rxn_name[r] <- rxns.lut$name[sel] )
        ( df.Rxns$rxn_eqn[r] <- rxns.lut$equation[sel] )
        ( df.Rxns$rxn_defn[r] <- rxns.lut$definition[sel] )
        
        # extract compound info
        
        #df.Rxns$rxn_eqn[r]
        #[1] "(1) cpd00010[0] + (1) cpd29672[0] <=> (1) cpd00045[0] + (1) cpd11493[0]"
        #[1] "(45) cpd00144[0] + (45) cpd00175[0] <=> (45) cpd00014[0] + (45) cpd00091[0] + (1) cpd15634[0]"
        
        ( compds.idx <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd")[[1]][,"start"] )
        # 5 23 43 61
        # 6 25 46 65 83
        
        ( compds <- as.character( str_extract_all(string = df.Rxns$rxn_eqn[r], pattern = "cpd.....", simplify = TRUE) ) )
        # "cpd00010" "cpd29672" "cpd00045" "cpd11493"
        
        if (length(compds)>=1) {
          
          df.Rxns$compds[r] <- paste0(compds, collapse = ";")
          
          ## get compound coefficients?
          start_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\(")[[1]][,"start"]
          end_brackets <- str_locate_all(string = df.Rxns$rxn_eqn[r], pattern = "\\)")[[1]][,"start"]
          ( compd.coeff <- as.numeric( substring(text = df.Rxns$rxn_eqn[r], first = start_brackets+1, last = end_brackets-1)) )
          
          df.Rxns$compd_coef[r] <- paste0(compd.coeff, collapse = ";")
          
          # get formulas of compounds
          
          formx <-filter(compounds.lut, id %in% compds )
          row.names(formx) <- formx$id
          ( formx.char <- formx[compds, ]$formula )
          # "C21H32N7O16P3S" "HOR"            "C10H11N5O10P2"  "C11H22N2O7PRS" 
          # "C15H19N2O18P2"      "C17H25N3O17P2"      "C9H12N2O12P2"       "C9H11N2O9P"         "C630H945N45O630P45"
          # "C7H7O7" "H2O"    "C7H5O6"
          df.Rxns$chem_formx[r] <- paste0(formx.char, collapse = ";")
          
          ( compd.names <- formx[compds, ]$name )
          # "2-methyl-trans-aconitate" "cis-2-Methylaconitate"
                    
          temp.df.Compounds <- data.frame(superfocus_fxn=fxn.superfocus.rowlabel,f=f, f__in=f.in,rxn_id= f.in.rxns[r], 
                                          cpd_id=compds, cpd_name=compd.names, cpd_form=formx.char, cpd_molar_prop=compd.coeff #, 
                                          #OC_x=OC_ratio, HC_y=HC_ratio , NC_z=NC_ratio 
          )
          
        } else {
          # No specified reaction equation or chemical formula info
          df.Rxns$compds[r] <- NA
          df.Rxns$compd_coef[r] <- NA
          df.Rxns$chem_formx[r] <- NA
          
          temp.df.Compounds <- NA
          
        }
        
        if (r==1) { df.Compounds <- temp.df.Compounds }
        
        if (r>1 & is.data.frame(df.Compounds) & is.data.frame(temp.df.Compounds)) { df.Compounds <- rbind(df.Compounds, temp.df.Compounds) }
        
        # clean up - if there are additional reactions?
        temp.df.Compounds <- NA
        
      } # END loop for r - rxn_id's per f/f.in
      
    } # END else loop when reactions identified
    
    # store results corresponding to each sub-reaction of each Superfocus row
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "fxns" ]] <- df.fxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "rxns" ]][[ f ]] <- df.Rxns
    
    if (f==1) { fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]] <- list() } # set this only once
    fxn.list[[ fxn.superfocus.rowlabel  ]][[ "compounds" ]][[ f ]] <- df.Compounds
    
    
  } # END loop - f in 1:length(fxns)) - to account for multiple functions/reactions reported in each row of Superfocus outputs
  
  saveRDS(object = fxn.list, file = paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds") )
    
} # END function to be run in parallel for each superfocus row


# # # # # # # # # # # # # # # # # #

no_forks <- 8

# this makes clusters on Unix-like system (may need to use other alternative for Windows)
cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
registerDoParallel(cl)

foreach(i=1:dim(df.tax)[1] , .packages=c('stringr', 'dplyr')) %dopar%  #
  get_rxns_and_compounds_indiv( df.tax=df.tax, subsys.lut=subsys.lut, rxns.lut=rxns.lut, rxn_pathways.lut=rxn_pathways.lut )

stopCluster(cl)


message("\n## assemble results")

message("\n(num_results_files <- dim(df.tax)[1])")
(num_results_files <- dim(df.tax)[1])

# assemble all compound data outputs
# start with blank row

df.out <- data.frame(superfocus_fxn=NA, f=NA, f__in=NA, rxn_id=NA, cpd_id=NA, cpd_name=NA, cpd_form=NA, cpd_molar_prop=NA )

for (i in 1:num_results_files) {
  fxn.superfocus.rowlabel <- row.names(df.tax)[i]
  temp <- readRDS(paste0(temp_dir,"/fxn-list-",fxn.superfocus.rowlabel,".rds"))
  
  f_no <- length( temp[[1]][["compounds"]] )
  
  for (f in 1:f_no) {
    #f<-2
    # only add non-NA results
    if (is.data.frame( temp[[1]][["compounds"]][[f]] )) {
      
      df.temp <- temp[[1]][["compounds"]][[f]]
      ok <- complete.cases(df.temp)
      df.temp <- df.temp[ which(ok==TRUE), ] # updated version will include some compounds with vK coordinates that are NA. vK coordinates are considered later
      df.out <- rbind(df.out,df.temp)
    }
  }
  print(paste0("added df ",i," of ",num_results_files ))
  
}


message("\nstr(df.out)")
str(df.out)


saveRDS(object = df.out, file = paste0("df.out--get_rxns_and_compounds_indiv-",this_study,".RDS"))

# remove NA first row
message("\nhead(df.out)")
head(df.out)

df.out <- df.out[-1, ]

message("\ndim(df.out)")
dim(df.out)


message("\n## normalise molar_prop to cpd_relabun so total of 1 per superfocus function")

df.out$cpd_molar_prop_norm <- NA

message("\nlength(unique(df.out$superfocus_fxn))")
length(unique(df.out$superfocus_fxn))

message("\nphy")
phy

message("\n% of functions represented - with compound information")
100*(length(unique(df.out$superfocus_fxn)) / ntaxa(phy))


fxns_found <- unique(df.out$superfocus_fxn)

for (k in 1:length(fxns_found)) {
  #k<-1
  this_fxn <- fxns_found[k]
  sel <- which(df.out$superfocus_fxn == this_fxn)
  
  sum_molar_prop <- sum( df.out$cpd_molar_prop[sel], na.rm = TRUE)
  # calculate 
  
  df.out$cpd_molar_prop_norm[sel] <- df.out$cpd_molar_prop[sel]/sum_molar_prop
  
  print(paste0("completed ",k))
  
}

message("\nsum(df.out$cpd_molar_prop_norm)")
sum(df.out$cpd_molar_prop_norm)

message("\nsample_sums(phy)")
sample_sums(phy)

message("\ngetwd()")
getwd()

saveRDS(object = df.out, file = paste0("df.out--tidy-compounds_indiv-cpp3d-",this_study,".RDS"))



message("\n### 2) get cpd rel abun per sample")
message("\n# # # # # # # # # #")


df.OTU <- as.data.frame( phy@otu_table ) # this is Superfocus functional relative abundance data represented in phyloseq OTU abundance table
message("\ndim(df.OTU)")
dim(df.OTU)


get_cpd_relabun_per_sample <- function(phy_in, dat.cpd) {
    
  this_samp <- sample_names(phy_in)[i]
  df.OTU <- as.data.frame( phy_in@otu_table[ ,this_samp] )
  
  dat.cpd$sample <- this_samp
  
  dat.cpd$cpd_rel_abun_norm <- NA
  
  fxns_all <- row.names(df.OTU)
  
  for (k in 1:length(fxns_all)) {
    #k<-1
    this_fxn <- fxns_all[k]
    sel <- which(dat.cpd$superfocus_fxn == this_fxn)
    
    if (length(sel)>=1) {
      dat.cpd$cpd_rel_abun_norm[sel] <- df.OTU[this_fxn, ]*dat.cpd$cpd_molar_prop_norm[sel]
      
    }
  } # END rel abun values for all relevant functions added
  
  saveRDS(object = dat.cpd, file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
  
} # END


no_forks <- 8

# this makes clusters on Unix-like system
cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
registerDoParallel(cl)

foreach(i=1: length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  get_cpd_relabun_per_sample( phy_in = phy, dat.cpd = df.out)

stopCluster(cl)


message("\n## assemble results")

# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
head(dat)

for ( i in 2:length(sample_names(phy)) ) {
  this_samp <- sample_names(phy)[i]
  temp <- readRDS( file = paste0(temp_dir,"/dat.cpd-",this_samp,".rds") )
  dat <- rbind(dat, temp)
  print(paste0("completed ",i))
}


saveRDS(object = dat, file = paste0("dat.cpd-long-all-samps-cpp3d-",this_study,".rds") )

rm(temp)

message("\nstr(dat)")
str(dat)

message("\nsum(dat$cpd_rel_abun_norm)")
sum(dat$cpd_rel_abun_norm)

message("\naverage functional relative abundance per sample")
message("\nsum(dat$cpd_rel_abun_norm)/nsamples(phy)")
sum(dat$cpd_rel_abun_norm)/nsamples(phy)

message("\nnames(dat)")
names(dat)

message("\nlength(unique(dat$cpd_id))")
length(unique(dat$cpd_id))




message("\n### 3) collate_compounds within each sample")
message("\n# # # # # # # # # #")


unique_cpd <- unique(dat$cpd_id)
samp_names <- sample_names(phy)


collate_compounds <- function(dat.cpd, unique_cpd, samp) {
  #i<-1
  #samp = samp_names[i]
  #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
  
  this_samp <- samp
  
  cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
  
  for (c in 1:length(unique_cpd)) {
    #c<-1
    this_cpd <- unique_cpd[c]
    sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    if (length(sel.cpd) >=1) {
      cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
    }
    
  } # END all compounds
  
  saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
  
} # END



no_forks <- 4

# this makes clusters on Unix-like system
cl<-makeForkCluster(nnodes = no_forks)   # no of nodes will depend on your HPC facility
registerDoParallel(cl)

foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])

stopCluster(cl)


message("\n## assemble results")

# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat.cpd.collate <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
head(dat.cpd.collate)

for ( i in 2:length(sample_names(phy)) ) {
  this_samp <- sample_names(phy)[i]
  temp <- readRDS( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") )
  
  dat.cpd.collate <- rbind(dat.cpd.collate, temp)
  
  print(paste0("completed ",i))
}


message("\nstr(dat.cpd.collate)")
str(dat.cpd.collate)

message("\nsum(dat.cpd.collate$cpd_rel_abun)")
sum(dat.cpd.collate$cpd_rel_abun)

message("\nsum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))")
sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))

saveRDS(object = dat.cpd.collate, file = paste0("dat.cpd.collate-all-samps-cpp3d-",this_study,".rds" ))

# END
