
# ======================================================================
# Set the random number seed as assigning individual ages with the ALK
# has a random component. This will increase reproducibility.
set.seed(9343888)



# ======================================================================
# Load the data.frames and do initial wranglings
## Load WBIC characteristics
wbic <- read.csv("data/prepped/wbicInfo.csv",stringsAsFactors=FALSE)
## Load length data
fmdb <- read.csv("data/prepped/fmdb_WAE.csv", stringsAsFactors = F) 
#fmdb = read.csv('results/age_lw/age_lw_all.csv', stringsAsFactors = F)
## Population estimates
### Removed WBIC_YEARs for which no FMDB data existed
pe <- read.csv("data/prepped/PE.csv")
rows2delete <- which(!pe$wbic_year %in% unique(fmdb$wbic_year))
pe <- pe[-rows2delete,]
cat(length(rows2delete),"WBIC_YEARs were removed from PE data.frame
    because no matching data in FMDB data.frame.")
## Load age-length-key information
### Only retain information for ALKs that are valid to use
### Retain only variables that are need when using the ALK
ALKInfo <- read.csv("data/prepped/ALKInfo.csv",stringsAsFactors=FALSE) %>%
  filterD(use=="yes") %>%
  select(type,which,ename,sname)
## Load weight-length regression results
### Only retain regressions results that are valid to use
### Remove variables that defined use and reason for not using
LWRegs <- read.csv("data/prepped/LWregs.csv",stringsAsFactors=FALSE) %>%
  filterD(use=="yes") %>%
  select(-use,-reason)



# ======================================================================
## Prepare for P/B calculating loop
wys <- as.character(unique(fmdb$wbic_year))
ttl.wys <- length(wys)
reg.src <- reg.type <- alk.src <- alk.type <- alk.note <- character(ttl.wys)
PE <- HA <- numAges <- n <- P <- B <- numeric(ttl.wys)
minAge <- maxAge <- numAgeGaps <- maxMissingAges <- numeric(ttl.wys)
minAgeR <- maxAgeR <- minLen <- numNegP <- numeric(ttl.wys)
maxMissingAgeAtEnd <- character(ttl.wys)
