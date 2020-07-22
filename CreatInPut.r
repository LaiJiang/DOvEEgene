
#####################################################################################
#an example to create data and make predictions based on
#existing codes and fitted model
# rm(list=ls())
#######################################define working directory

PATH.work <- "~/Ragousis/pipeline/NewVariants/"
#PATH.work <- "~/Documents/doveegene/CodeNov2019/pipeline_pred_modified/"
PATH.dat <- paste0(PATH.work, "dat/")

PATH.func <- paste0(PATH.work, "func/")

PATH.save <- paste0(PATH.work, "save/")

#######################################load raw data

#load clinical information of mutation carriers
clinc.info <- read.csv( paste0(PATH.dat, "new_clinical.csv"))

#load patients-mutations links
sample.mut.link <-  read.csv( paste0(PATH.dat, "ind_data_new_filter.csv"))

#the VAF value=1 will cause trouble in log(1-VAF). Thus we set 1 ~= 0.999.
#note the maximum VAF value in previous filter = 0.95. 
sample.mut.link$VAF[sample.mut.link$VAF==1] <- 0.999
  
#load mutation information/annotations
mut.info <- read.csv(paste0(PATH.dat, 'var_data_new_filter.csv'))

##########################################clear data
#check factor variables and convert to characters (since we need to count the number of indels..etc)
factor.var <- sapply(mut.info, is.factor)
if(length(factor.var)>0)mut.info[factor.var] <- lapply(mut.info[factor.var], as.character)

#missing values N/A replaced by 0
#as long as 0 values are consistenly on the other side (< or >) of all non-missing cases 
#then this coding would have no impact on analysis
mut.info[mut.info=="N/A"] <- 0

mut.info$type[mut.info$type=="SNV"] <- "snp"
mut.info$type[mut.info$type!="snp"] <- "indel"

mut.info$impact_severity[mut.info$impact_severity=="MODERATE"] <- "MED"



#######################################load all functions for collapsing
source(paste0(PATH.func, "func_collapse.r"))


#List all functions for collapsing: name start with "fn."
all.func.names <- lsf.str()
all.func.names <- unique(all.func.names[grep("fn.",all.func.names)])

#collapse variants in var.data into ind.data using specified function list (for collapsing)
options(warn=0)

#create summarized data
t.all <- perpersonsummary(sample.mut.link, mut.info, 
                          base::mget(all.func.names) )

#patient IDs assign to each row
rownames(t.all)<-names(table(sample.mut.link$patientID))

#function names assign to each column 
colnames(t.all) <- all.func.names

#convert to dataframe
test.t.all <- as.data.frame(t.all)

#save the collaped data
save(test.t.all, file=paste0(PATH.save, "all_collapse_2207.RData"))


  