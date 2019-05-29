## Clear workspace and console
rm(list=ls()); cat("\014")

# ######################################################################
# ======================================================================
# This script computes weight-length regressions by WBIC_YEAR, WBIC,
# lake classification, and overall (no grouping) and writes a file with
# the regression coefficients and other regression information. This 
# script need not be run again unless the underlying weight-length data
# are modified. The file created will be used to predict weights from
# lengths for fish that were not weighed.
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
library(nlme)
source("code/helpers/productionHelpers.R")

# ======================================================================
# Load length-weight data
##  Added logged versions of length and weight
lw <- read.csv("data/prepped/len_wt.csv") %>%
  mutate(logl=log(len.mm),logw=log(wt))


# ======================================================================
# Calculate length-weight regressions
## For each WBIC_YEAR
##   na.action needed because of NAs in other variables in data.frame
lw_wy <- nlme::lmList(logw~logl|wbic_year,data=lw,na.action=na.pass)
wy_res <- LWINFO(lw_wy,"WBIC_YEAR")
## Try this to see an individual regression plot (can change number)
fitPlot(lw_wy[[names(lw_wy)[1]]],main=names(lw_wy)[1])

## For each WBIC
lw_w <- nlme::lmList(logw~logl|wbic,data=lw,na.action=na.pass)
w_res <- LWINFO(lw_w,"WBIC")
fitPlot(lw_w[[names(lw_w)[1]]],main=names(lw_w)[1])

## For each lake class
lw_c <- nlme::lmList(logw~logl|class,data=lw,na.action=na.pass)
c_res <- LWINFO(lw_c,"CLASS")
fitPlot(lw_c[[names(lw_c)[1]]],main=names(lw_c)[1])

## For all fish
lw_a <- lm(logw~logl,data=lw)
tmp <- summary(lw_a)
tmp2 <- model.frame(lw_a)$logl
a_res <- data.frame(coef(lw_a)[1],coef(lw_a)[2],rsq=tmp$r.squared,
                    n=tmp$df[2]+2,minLen=exp(min(tmp2)),
                    maxLen=exp(max(tmp2)),which="ALL",type="ALL") %>%
  rename(loga=coef.lw_a..1.,b=coef.lw_a..2.)
fitPlot(lw_a)


# ======================================================================
# Combine all regressions into one data.frame
## Add a "use" variable
##    NO for "small" n, "low" r-squared, weird b ("far" from 3)
##    as defined by these "cut"offs
n.cut <- 25; rsq.cut <- 0.85; b.cut <- c(2,4)
## Add exponentiated intercept (a) variable
lw_res <- rbind(wy_res,w_res,c_res,a_res) %>%
  mutate(use=case_when(n < n.cut ~ "NO",rsq < rsq.cut ~ "NO",
                       b < b.cut[1] | b > b.cut[2] ~ "NO",
                       TRUE ~ "yes"),
         reason=case_when(n < n.cut ~ paste("n<",n.cut),
                          rsq < rsq.cut ~ paste("r^2<",rsq.cut),
                          b < b.cut[1] | b > b.cut[2] ~ "b far from 3",
                          TRUE ~ ""),
         use=factor(use),reason=factor(reason),
         a=exp(loga)) %>%
  select(use,type,which,loga,a,b,n,rsq,minLen,maxLen,reason)
headtail(lw_res)


# ======================================================================
# Output results to a file in data/prepped/
if (TRUE) write.csv(lw_res,"data/prepped/LWregs.csv",
                                 quote=FALSE,row.names=FALSE)
