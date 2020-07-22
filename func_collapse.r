##############################################################################
# This code only saves functions to collapse variables by patient
# This file will be called in creat_dat.r
##############################################################################

#summary by person using a specified function list: functioncalls
perpersonsummary <- function(fn.ind.data, fn.var.data, functioncalls) {

    patientids <- table(fn.ind.data$patientID)
  for (ii in (1:length(patientids))) {
    one.data <- fn.ind.data[fn.ind.data$patientID==names(patientids)[ii],]
    variant.data <- fn.var.data[match(one.data$variantID, fn.var.data$variantID),]
    
    sone.data <- NULL
    for(func_ii in 1:length(functioncalls)){ 
      sone.data <- c(sone.data, functioncalls[[func_ii]](one.data, variant.data))  
    }
    
    if (ii==1) {summar.data <- sone.data}
    if (ii>1) {summar.data <- rbind(summar.data, sone.data)}
  }
  return(summar.data)
}




##############################################################################
#Below shows all functions that will be called into perpersonsummary() by creat_dat.r

fn.mnvaf <- function(one.data, variant.data) {return(mean(one.data$VAF)) }
fn.maxvaf <- function(one.data, variant.data) {return(max(one.data$VAF))}

fn.max.polyphen <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$polyphen_score), na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
  }

fn.max.sift <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$sift_score), na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
  }

fn.max.gerp_bp <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$gerp_bp_score), na.rm=TRUE)
  return(tmp)}

fn.max.cadd <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$cadd_scaled), na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
  }


fn.max.clinPred <- function(one.data, variant.data)  {
  tmp <- max(as.numeric(variant.data$clinPred), na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
}


fn.mn.polyphen <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score)*one.data$VAF,na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
  
  }

fn.1me.polyphen <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score) *
               log(1-one.data$VAF), na.rm=TRUE)
  return( ifelse(  is.na(tmp), 0, tmp  ) )
  }



fn.maxvaf.polyphen <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * one.data$VAF
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)}


fn.mn.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }

fn.1me.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score) *
               log(1-one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}



fn.maxvaf.sift <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$sift_score) * one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp)}

fn.mn.gerp_bp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }

fn.1me.gerp_bp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score) *
               log(1-one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}

fn.maxvaf.gerp_bp <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$gerp_bp_score) * one.data$VAF,na.rm=TRUE)
  return(tmp)}

fn.mn.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }

fn.1me.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled) *
               log(1-one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}



fn.maxvaf.cadd <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$cadd_scaled)*one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp)}


fn.mn.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred)*one.data$VAF,na.rm=TRUE)
  if (is.na(tmp)) tmp <-0; return(tmp) }

fn.1me.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred) *
               log(1-one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}


fn.maxvaf.clinPred <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$clinPred) * one.data$VAF, na.rm=TRUE)
  if (is.na(tmp)) tmp <-0
  return(tmp) }


fn.number <- function(one.data, variant.data)  {return(nrow(one.data))}

fn.number.high <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='HIGH',])) }

fn.number.med <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='MED',])) }

fn.number.low <- function(one.data, variant.data) {
  return(nrow(one.data[variant.data$impact_severity=='LOW',])) }


fn.miss <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$type=='indel']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(length(tmp))}
}

fn.count.indel <- function(one.data, variant.data) {
  tmp<-sum(variant.data$type=="indel") 
  if(is.na(tmp))tmp<-0
  return(tmp)
}

####################################
#B count #SNPs
fn.count.snp <- function(one.data, variant.data) {
  f<-sum(variant.data$type=="snp") 
  if(is.na(f))f<-0
  return(f)
}

####################################
#C count total #alterrations
fn.count.all <- function(one.data, variant.data) {
  f<-sum(variant.data$type %in% c("snp" ,"indel" ) ) 
  if(is.na(f))f<-0
  return(f)
}



####################################
#VAF weighted A
fn.VAF.count.indel <- function(one.data, variant.data) {
  f<-sum( (variant.data$type=="indel") * one.data$VAF) 
  if(is.na(f))f<-0
  return(f)
}



####################################
#VAF weighted B
fn.VAF.count.snp <- function(one.data, variant.data) {
  f<-sum( (variant.data$type=="snp")* one.data$VAF )
  if(is.na(f))f<-0
  return((f))
}



####################################
#VAF weighted c
fn.VAF.count.all <- function(one.data, variant.data) {
  f<-sum( (variant.data$type %in% c("snp" ,"indel" ) ) * one.data$VAF )
  if(is.na(f))f<-0
  return(f)
}

####################################
#maximum (gerp * vaf *A)

fn.max.vaf.gerp.indel <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="indel") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

####################################
#maximum (gerp * vaf *c)

fn.max.vaf.gerp.all <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type%in%c("indel","snp") ) 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

####################################
#mean (gerp * vaf *A)

fn.mn.vaf.gerp.indel <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="indel") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

####################################
#mean (gerp * vaf *c)

fn.mn.vaf.gerp.all <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type%in%c("indel","snp") ) 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

####################################

#maximum (polyphen score * vaf * B)
fn.max.vaf.poly.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


####################################
fn.max.vaf.sift.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$sift_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


#I have minus infinity here.

####################################
fn.max.vaf.gerp.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}
#########
fn.max.vaf.cadd.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$cadd_scaled) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


####################################
fn.max.vaf.clin.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$clinPred) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- max(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}

####################################

fn.mean.vaf.poly.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$polyphen_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


#I have minus infinity here.
####################################
fn.mean.vaf.sift.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$sift_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}
############################
fn.mean.vaf.gerp.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$gerp_bp_score) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


####################################
fn.mean.vaf.cadd.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$cadd_scaled) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


####################################
fn.mean.vaf.clin.snp <- function(one.data, variant.data) {
  tmp <- as.numeric(variant.data$clinPred) * 
    one.data$VAF * (variant.data$type=="snp") 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}


####################################

fn.no.variant <- function(one.data, variant.data) {
  tmp <- as.numeric(sum(!is.na(one.data$variantID)) )
  return(tmp)  
}




fn.VAF.high <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='HIGH']
  return( ifelse(  length(tmp)>0, sum(log(tmp)), 0  ) )
}

fn.VAF.med <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='MED']
  return( ifelse(  length(tmp)>0, sum(log(tmp)), 0  ) )
}

fn.VAF.low <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='LOW']
  return( ifelse(  length(tmp)>0, sum(log(tmp)), 0  ) )
}

fn.log.OM.VAF <- function(one.data, variant.data) {return(sum(log(1-one.data$VAF)))  }

fn.entropy <-  function(one.data, variant.data) {
  return(sum(one.data$VAF*log(one.data$VAF)))  }


fn.E.high <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='HIGH']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
  if(is.na(tmp)==1){return(0)}
}

fn.E.med <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='MED']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
}


fn.E.low <- function(one.data, variant.data) {
  tmp <- one.data$VAF[variant.data$impact_severity=='LOW']
  if (length(tmp)==0) {return(0)}
  if (length(tmp)>0) {return(sum(tmp*log(tmp)))}
}


fn.OME <- function(one.data, variant.data) {
  return(sum((1-one.data$VAF)*log(1-one.data$VAF)))  }


fn.max.size.base <- function(one.data, variant.data) {
  tmp <- max(as.numeric(variant.data$end - variant.data$start), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}

fn.sum.size.base <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$end - variant.data$start), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}



fn.max.size.bp.change <- function(one.data, variant.data) {
  tmp <- max(-( nchar(variant.data$alt) - nchar(variant.data$ref)), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}



fn.sum.size.bp.change <- function(one.data, variant.data) {
  tmp <- sum(-( nchar(variant.data$alt) - nchar(variant.data$ref)), na.rm=TRUE)
  #if (is.na(tmp)) tmp<--20  # there are no missing values
  return(tmp)}


fn.mean.VAF.size.change  <- function(one.data, variant.data) {
  tmp <- as.numeric( ( nchar(variant.data$alt) - nchar(variant.data$ref)) ) * 
    one.data$VAF 
  tmp2 <- mean(tmp, na.rm=TRUE)
  if (is.na(tmp2)) tmp2 <- 0; return(tmp2)
  
}





fn.sum.log.vaf <- function(one.data, variant.data) {
  return(sum( log(one.data$VAF)) ) }

fn.sum.log.om.vaf <- function(one.data, variant.data) {
  return(sum( log(1-one.data$VAF)) ) }


fn.entropy.cadd <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$cadd_scaled) *
               log(one.data$VAF)*(one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}



fn.entropy.clinPred <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$clinPred) *
               log(one.data$VAF)*(one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}


fn.entropy.poly <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$polyphen_score) *
               log(one.data$VAF)*(one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}


fn.entropy.sift <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$sift_score) *
               log(one.data$VAF)*(one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}

fn.entropy.gerp <- function(one.data, variant.data) {
  tmp <- sum(as.numeric(variant.data$gerp_bp_score) *
               log(one.data$VAF)*(one.data$VAF), na.rm=TRUE)
    return( ifelse(  is.na(tmp), 0, tmp  ) )}

fn.cosm.id <- function(one.data, variant.data) {
  return( as.numeric( nrow(one.data) - nrow(one.data[variant.data$cosmic_ids=='None',])>0) ) }


fn.age <-  function(one.data, variant.data) {
  return( as.numeric( mean(one.data$age) ) ) }
