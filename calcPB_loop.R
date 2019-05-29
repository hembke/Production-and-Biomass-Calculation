## Clear workspace and console
rm(list=ls()); cat("\014")
# Load required packages
library(FSA); library(dplyr); library(magrittr);library(zoo)
# Source local helper files
source("code/helpers/calcPB.R")
source("code/helpers/productionHelpers.R")

# ######################################################################
# ======================================================================
# This script will loop through all valid WBIC_YEARs (i.e., has a PE and
# has some measured fish in the fmdb) and calculate P and B for each.
# Those values, along with some other information are then output to a
# file in the results/ folder for further analysis.
#
# This script requires that the Data_Prepper, calcLWRegs, and calcALKs
# scripts have all been successfully run (i.e., their resultant files
# were created and stored in the data/prepped/ folder.)
# ======================================================================
# ######################################################################

# ======================================================================
# Analysis choices that can be made (at this stage)
## Minimum size of sample when computing P & B
n.cut <- 30
## Minimum number of ages in sample when computing P & B
ages.cut <- 5
## Type of ALK to use ("empirical" or "smoothed")
alk2use <- "smoothed"
## Assign individual weights (TRUE) or est. mean weight from mean length
useIndivWeights <- TRUE
## An optional name to add as a suffix to the results file. Use to keep
## results separate based on choices made above. Use NULL for no suffix.
results.suffix <- "SM_IW_2017"
## Summary plots of PB calculations can be made, but this slows the loop
## considerably. Perhaps turn this off until near final run.
makeSumPlots <- FALSE


# ======================================================================
# Load CalcPB_Setup.R script
source("code/helpers/calcPB_Setup.R")


# ======================================================================
# Compute P & B (et al.) for each WBIC_YEAR
if (useIndivWeights) {
  for (i in 1:ttl.wys) {
    print(paste0("i=",i," of ",ttl.wys," (",wys[i],")"))
    fmdb_1 <- filterD(fmdb,wbic_year==wys[i])
    # Get WBIC_YEAR PE and WBIC size
    PE[i] <- getPE(fmdb_1,pe)
    HA[i] <- getSize(fmdb_1,wbic)
    # Apply ALKs
    tmp <- doALK(fmdb_1,ALKInfo,alk2use)
    fmdb_1 <- tmp$df
    alk.src[i] <- tmp$which; alk.type[i] <- tmp$type; 
    # Apply LWRegs (must check fmdb_1 as it may have no fish after ALK)
    if (nrow(fmdb_1)>0) {
      tmp <- doLWReg(fmdb_1,"len.mm","wt",LWRegs)
      fmdb_1 <- tmp$df
      reg.src[i] <- tmp$reg$which; reg.type[i] <- tmp$reg$type    
    }
    # Get overall sample size
    n[i] <- nrow(fmdb_1)
    # Get number of ages present
    numAges[i] <- length(unique(fmdb_1$age))
    # Get minimum length
    minLen[i] <- min(fmdb_1$len.mm)
    # Calculate P and B (if more than one age-class)
    ## get number and mean weight (kg) in each age-class
    ## expand number in sample to number in popn (with PE value)
    ## add a total weight
    ## send to calcPB to compute P and B
    if (numAges[i]>1) {
      sum_1 <- group_by(fmdb_1,age) %>%
        summarize(snum=n(),mwt=mean(wt)/1000) %>%
        mutate(pnum=snum/sum(snum)*PE[i],twt=pnum*mwt) %>%
        calcPB(age.c="age",num.c="pnum",twt.c="twt",
               area=HA[i],adjAgeGaps=TRUE,lbl=wys[i])
      P[i] <- sum_1$P; B[i] <- sum_1$B
      # get some characteristics of ages present
      minAge[i] <- min(sum_1$df$age); maxAge[i] <- max(sum_1$df$age)
      ## Number of gaps in ages (e.g., 8 and 10, but not no 9)
      tmp <- diff(sum_1$df$age)
      numAgeGaps[i] <- sum(tmp>1)
      ## How many ages are missing in the largest gap?
      maxMissingAges[i] <- max(tmp)-1
      ## Is the largest gap near the oldest ages (within 3 age-classes)?
      ## If largest gap occurs more than once, then say NOT at end
      tmp2 <- which(tmp==(maxMissingAges[i]+1))
      if (length(tmp2)==1) maxMissingAgeAtEnd[i] <- 
        ifelse(tmp2 %in% (length(tmp)-2):length(tmp),"yes","")
      ## Describe longest run of continuous ages
      tmp <- split(sum_1$df$age,cumsum(c(TRUE,diff(sum_1$df$age)!=1)))
      tmp <- tmp[[which.max(sapply(tmp,length))]]
      minAgeR[i] <- min(tmp); maxAgeR[i] <- max(tmp)
      # How many age-class Ps were negative
      numNegP[i] <- sum(sum_1$df$P<0,na.rm=TRUE)
      ## Write out the calculation table so it can be examined later
      write.csv(sum_1$df,paste0("results/CalcPB_Tables/PB",
                                ifelse(is.null(results.suffix),"","_"),
                                results.suffix,
                                ifelse(is.null(results.suffix),"","_"),
                                wys[i],".csv"),quote=FALSE,row.names=FALSE)
      ## Also put out a summary graphic (if asked for ... very slow)
      if (makeSumPlots) {
        jpeg(paste0("results/CalcPB_Graphs/PB",
                    ifelse(is.null(results.suffix),"","_"),
                    results.suffix,
                    ifelse(is.null(results.suffix),"","_"),
                    wys[i],".jpg"),width=3*480,res=144)
        plot(sum_1)
        dev.off() 
      }
    } else P[i] <- B[i] <- NA
  }
} else { ## compute mean weights from mean lengths-at-age
  for (i in 1:ttl.wys) {
    print(paste0("i=",i," of ",ttl.wys," (",wys[i],")"))
    fmdb_1 <- filterD(fmdb,wbic_year==wys[i])
    # Get WBIC_YEAR PE and WBIC size
    PE[i] <- getPE(fmdb_1,pe)
    HA[i] <- getSize(fmdb_1,wbic)
    # Apply ALKs
    tmp <- doALK(fmdb_1,ALKInfo,alk2use)
    fmdb_1 <- tmp$df
    alk.src[i] <- tmp$which; alk.type[i] <- tmp$type; alk.note[i] <- tmp$note
    # Get overall sample size
    n[i] <- nrow(fmdb_1)
    # Get number of ages present
    numAges[i] <- length(unique(fmdb_1$age))
    # Get minimum length
    minLen[i] <- min(fmdb_1$len.mm)
    # Calculate P and B (if more than one age-class)
    ## get number and mean, SD and CV weight (g) in each age-class
    ## correct the mean weight using eqn 16 from Beyer (1991)
    ## convert mean weights to kg
    ## expand number in sample to number in popn (with PE value)
    ## add a total weight
    ## send to calcPB to compute P and B
    if (numAges[i]>1) {
      sum_1 <- group_by(fmdb_1,age,wbic_year,wbic,year,class) %>%
        summarize(snum=n(),mlen=mean(len.mm),sdlen=sd(len.mm)) %>%
        mutate(cvlen=sdlen/mlen) %>%
        ungroup() %>% as.data.frame()
      tmp <- doLWReg(sum_1,"mlen","mwt",LWRegs)
      b <- LWRegs$reg$b     # isolate b coefficient
      sum_1 <- tmp$df %>%
        select(-wbic_year,-wbic,-year,-class) %>%
        mutate(mwt=mwt*((1+cvlen^2)^(b*(b-1)/2)),
               mwt=mwt/1000) %>%
        mutate(pnum=snum/sum(snum)*PE[i],twt=pnum*mwt) %>%
        calcPB(age.c="age",num.c="pnum",twt.c="twt",area=HA[i],
               adjAgeGaps=TRUE,lbl=wys[i])
      P[i] <- sum_1$P; B[i] <- sum_1$B
      # get characteristics of LW used
      reg.src[i] <- tmp$reg$which; reg.type[i] <- tmp$reg$type
      # get some characteristics of ages present
      minAge[i] <- min(sum_1$df$age); maxAge[i] <- max(sum_1$df$age)
      ## Number of gaps in ages (e.g., 8 and 10, but not no 9)
      tmp <- diff(sum_1$df$age)
      numAgeGaps[i] <- sum(tmp>1)
      ## How many ages are missing in the largest gap?
      maxMissingAges[i] <- max(tmp)-1
      ## Is the largest gap near the oldest ages (within 3 age-classes)?
      ## If largest gap occurs more than once, then say NOT at end
      tmp2 <- which(tmp==(maxMissingAges[i]+1))
      if (length(tmp2)==1) maxMissingAgeAtEnd[i] <- 
        ifelse(tmp2 %in% (length(tmp)-2):length(tmp),"yes","")
      ## Describe longest run of continuous ages
      tmp <- split(sum_1$df$age,cumsum(c(TRUE,diff(sum_1$df$age)!=1)))
      tmp <- tmp[[which.max(sapply(tmp,length))]]
      minAgeR[i] <- min(tmp); maxAgeR[i] <- max(tmp)
      ## How many age-class Ps were negative
      numNegP[i] <- sum(sum_1$df$P<0,na.rm=TRUE)
      ## Write out the calculation table so it can be examined later
      write.csv(sum_1$df,paste0("results/CalcPB_Tables/PB",
                                ifelse(is.null(results.suffix),"","_"),
                                results.suffix,
                                ifelse(is.null(results.suffix),"","_"),
                                wys[i],".csv"),quote=FALSE,row.names=FALSE)
    } else P[i] <- B[i] <- NA
  }
}


# ======================================================================
# Put Results together
split.wy <- do.call(rbind,strsplit(wys,"_"))
PB_res <- data.frame(wbic_year=wys,wbic=split.wy[,1],year=split.wy[,2],
                     n,minLen,numAges,minAge,maxAge,
                     numAgeGaps,maxMissingAges,maxMissingAgeAtEnd,
                     minAgeR,maxAgeR,
                     PE,HA,P,B,numNegP,
                     reg.type,reg.src,alk.type,alk.src,alk.note)
## Add use and reason variables
PB_res %<>% mutate(use=case_when(n<n.cut ~ "NO",
                                 numAges<ages.cut ~ "NO",
                                 (numAgeGaps>=3 & (maxAgeR-minAgeR<4)) ~ "NO",
                                 TRUE ~ "yes"),
                   reason=case_when(n<n.cut ~ paste0("n<",n.cut),
                                    numAges<ages.cut ~ paste0("numAges<",ages.cut),
                                    (numAgeGaps>=3 & (maxAgeR-minAgeR<4)) ~ "Age gaps issue",
                                    TRUE ~ ""))


# ======================================================================
## Write the file out to the results folder with a time stamp in name
write.csv(PB_res,paste0("results/PB",ifelse(is.null(results.suffix),"","_"),
                        results.suffix,".csv"),quote=FALSE,row.names=FALSE)
