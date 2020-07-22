##############
#example script for making  predictions

rm(list=ls())
#######################################define working directory
PATH.work <- "/mnt/GREENWOOD_JBOD1/GREENWOOD_BACKUP/home/lai.jiang/Ragousis/pipeline/"

PATH.dat <- paste0(PATH.work,  "NewVariants/dat/")

PATH.func <- paste0(PATH.work, "NewVariants/func/")

PATH.save <- paste0(PATH.work, "NewVariants/save/")

#######################################load collaped data
#clinc.info <- read.csv(paste0(PATH.dat,"/Clini_mutpatients.csv"))

#load("save/all_collapse_08072020RData")
load(paste0(PATH.save, "all_collapse_2207.RData" ))
#######################################load trained model and trained dataset
load(paste0(PATH.save, "PhaseII_saved_model.RData" ))

clinc.info.new <- read.csv(paste0(PATH.dat,"/new_clinical.csv"))

######################


#creat test.x and test.y (if any) for prediction
#we predict on all cases in test.t.all
#test.x <- test.t.all
#test.y <- (clinc.info$Dx!="Benign")[match(rownames(test.t.all),clinc.info$patientID)]
#test.y <- as.factor(test.y)
#levels(test.y)  <- c("Benign", "Cancer")


#creat test.x and test.y (if any) for prediction
#we predict on all cases in test.t.all
test.x <- test.t.all
test.y <- (clinc.info.new$dx!="Benign")[match(rownames(test.t.all),clinc.info.new$patientID)]
test.y <- as.factor(test.y)
levels(test.y)  <- c("Benign", "Cancer")


#######standradize test.x according to mean and sd 
#note that the standardization here is based on test.x alone
#Moreover, we can standarize test.x according to train.x (from which we trained model) if  
#both train.x  and test.x are obtained from the same platform (e.g. germline)

#Note: the standardization mn.x and sd.x should be loaded from the trained model
#because the test.x should be conformed to the distribution of train.x for staitstical validity.
sd.ts.x <- sweep(test.x, 2, mn.x, FUN = "-")
sd.ts.x <- sweep(sd.ts.x, 2, sd.x, FUN = "/")
###################################################
#now we predict test.x using the model trained before ()

pred.test <- predict(wsrf.spec.weights, newdata=sd.ts.x, type = "prob")

#plot predicted probability according to cancer/benign
plot(pred.test$Cancer,col=(test.y=="Cancer")+1, xlab="patients", ylab="Cancer probs",
     main="red:cancer")
abline(h=0.8)


#ROC curves

library(PRROC)

ROC_obj <- roc.curve(scores.class0 = pred.test$Cancer, 
                     weights.class0=(test.y=="Cancer"),
                       curve=TRUE)
plot(ROC_obj)

PR_obj <- pr.curve(scores.class0 = pred.test$Cancer,
                   weights.class0=(test.y=="Cancer"),
                   curve=TRUE)
plot(PR_obj)

