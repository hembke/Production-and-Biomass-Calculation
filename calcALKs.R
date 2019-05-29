## Clear workspace and console
rm(list=ls()); cat("\014")

# ######################################################################
# ======================================================================
# This script computes age-length keys by WBIC_YEAR, WBIC, lake
# classification, and overall (no grouping) and writes XXX. This 
# script need not be run again unless the underlying age-length data
# are modified. The files created will be used to predict ages from
# lengths for fish that were not aged.
#
# THIS SCRIPT SHOULD NOT NEED TO BE RUN AGAIN AS IT WILL OVERWRITE
# THE PREPPED FILE IN 'data/prepped/'.

# However, as a safeguard, the file will only be overwritten if the
# writePreppedFiles object below is set to TRUE.
writePreppedFiles <- TRUE
#
# This script requires that the Data_Prepper script has been
# successfully run (i.e., its resultant files were created and stored
# in the data/prepped/ folder.)
# ======================================================================
# ######################################################################

# ======================================================================
# Setup
## Load required packages
library(FSA)
library(dplyr)
library(magrittr)
library(nnet)
## A new function to construct ALKs and return some information from
## within a loop
ALKINFO <- function(ladat,group,type,writeFiles) {
  # Create raw two-way contingency table (precursor to empirical ALK)
  tmp <- xtabs(~lcat+age,data=ladat)
  # Get min, max, and width (i.e., should be minimum width among length
  # categories) of length categories
  tmp.lcats <- as.numeric(rownames(tmp))
  minLcat <- min(tmp.lcats)
  maxLcat <- max(tmp.lcats)
  w <- ifelse(length(tmp.lcats)>1,min(diff(tmp.lcats)),NA)
  # Get overall sample size
  n <- sum(tmp)
  # Create empirical ALK ... save out to R object
  alk <- prop.table(tmp,margin=1)
  fname1 <- paste("eALK",group,sep="_")
  if (writeFiles) save(alk,file=paste0("data/prepped/ALKs/",fname1,".RData"))
  # Create smoothed ALK ... save out to R object
  #   Restricted to those with numAges>3 and n>10 (small numbers of ages
  #   and small n were causing problems ... convergence but weird ALK)
  numAges <- ncol(tmp)
  if (numAges>3 & n>10) {
    # fit multinomial model ... capture.output used to hide iteration notes
    junk <- capture.output( mlr <- nnet::multinom(age~lcat,data=ladat,maxit=500) )
    # Check convergence ... 0 means converged, 1 means did not converge
    # If converged find smoothed ALK, if not do not
    if (mlr$convergence==0) {
      # create length categories over which to make the smootherd ALK
      # will create from minimum to maximum plus one category width. The
      # extra category width helps with fish larger than max observed len
      lens <- seq(minLcat,maxLcat+w,w)
      # Predicted probabilities of ages for each length category
      # I rounded to four decimal places to eliminate chances of a very
      # weird age prediction (i.e., without this there is a small chance
      # that, for example, a 200 mm fish could be predicted to be age 20,
      # or something similar).
      alk <- round(predict(mlr,data.frame(lcat=lens),type="probs"),3)
      # However, there were some problems with rows not summing to 1 
      # (e.g., summing to 1.01); thus force all rows to sum to 1.
      alk <- prop.table(alk,margin=1)
      # Give the row names the names of the lengths.
      rownames(alk) <- lens
      # Save the smoothed ALK to an R object
      fname2 <- paste("sALK",group,sep="_")
      if (writeFiles) save(alk,file=paste0("data/prepped/ALKs/",fname2,".RData"))      
    } else fname2 <- "NA"
  } else fname2 <- "NA"
  # return a 1 row data.frame with a bunch of information in it
  data.frame(which=group,type=type,n=n,
             minLen=minLcat,maxLen=maxLcat,numLens=nrow(tmp),
             minAge=min(as.numeric(colnames(tmp))),
             maxAge=max(as.numeric(colnames(tmp))),
             numAges=numAges,ename=fname1,sname=fname2,
             stringsAsFactors=FALSE)
}



# ======================================================================
# Load lengh-age data
## Added length categories
la <- read.csv("data/prepped/len_age.csv",stringsAsFactors=FALSE) %>%
  mutate(lcat=lencat(len.mm,w=20))



# ======================================================================
# Compute all ALKs and save out to individual R objects
## By WBIC_YEAR
wys <- unique(la$wbic_year)
res_wy <- ALKINFO(filterD(la,wbic_year==wys[1]),wys[1],
                  "WBIC_YEAR",writePreppedFiles)
for (i in 2:length(wys)) {
  print(paste0("i=",i," of ",length(wys)," for WBIC_YEARs"))
  res_wy <- rbind(res_wy,ALKINFO(filterD(la,wbic_year==wys[i]),
                                 wys[i],"WBIC_YEAR",writePreppedFiles))
}

## By WBIC
ws <- unique(la$wbic)
res_w <- ALKINFO(filterD(la,wbic==ws[1]),ws[1],"WBIC",writePreppedFiles)
for (i in 2:length(ws)) {
  print(paste0("i=",i," of ",length(ws)," for WBICs"))
  res_w <- rbind(res_w,ALKINFO(filterD(la,wbic==ws[i]),ws[i],"WBIC",writePreppedFiles))
}

## By lake class (must be patient)
lc <- unique(la$class)
lc <- lc[!is.na(lc)]
res_c <- ALKINFO(filterD(la,class==lc[1]),lc[1],"CLASS",writePreppedFiles)
for (i in 2:length(lc)) {
  print(paste0("i=",i," of ",length(lc)," for lake classes"))
  res_c <- rbind(res_c,ALKINFO(filterD(la,class==lc[i]),lc[i],
                               "CLASS",writePreppedFiles))
}

## All (must be very patient)
res_a <- ALKINFO(la,"ALL","ALL",writePreppedFiles)



# ======================================================================
# Combine all ALK information into one data.frame
## Add a "use" variable
##    NO for "small" n, few age-classes, few length-classes (as defined
##    by the "cut"offs below) or the multinomial model was not fit (i.e,
##    this would imply that the epirical ALK was questionable as well.)
n.cut <- 30; ages.cut <- 5; lens.cut <- 5
alk_res <- rbind(res_wy,res_w,res_c,res_a) %>%
  mutate(use=case_when(n < n.cut ~ "NO",numAges < ages.cut ~ "NO",
                       numLens < lens.cut ~ "NO",sname=="NA" ~ "NO",
                       TRUE ~ "yes"),
         reason=case_when(n < n.cut ~ paste("n <",n.cut),
                          numAges < ages.cut ~ paste("# ages <",ages.cut),
                          numLens < lens.cut ~ paste("# lcats <",lens.cut),
                          sname=="NA" ~ "Bad smoothed ALK",
                          TRUE ~ ""),
         use=factor(use),reason=factor(reason)) %>%
  select(use,type,which,n,minLen,maxLen,numLens,
         minAge,maxAge,numAges,ename,sname,reason)
headtail(alk_res)


# ======================================================================
# Output info results to a file in data/prepped/
if (writePreppedFiles) write.csv(alk_res,"data/prepped/ALKInfo.csv",
                                 quote=FALSE,row.names=FALSE)
