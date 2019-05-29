#' Calculate production and biomass with instantaneous growth method
#' 
#' Calculate production and biomass with instantaneous growth method when given a data.frame that contains abundance and total biomass for each age-class.
#' 
#' @param df A data.frame that contains variables that identify age-classes and abundance and total biomass for each age-class.
#' @param age.c A single numeric or string for the column in \code{df} that contains the age-classes.
#' @param num.c A single numeric or string for the column in \code{df} that contains the abundances for each age-classe.
#' @param age.c A single numeric or string for the column in \code{df} that contains the biomass for each age-class.
#' @param adjAgeGaps A logical that tells whether to adjust for missing age-classes (e.g, if some age-classes represent more than one age-class; see the Escanaba Lake example). When \code{TRUE} the biomass is divided by number of years between ages to put in a per year basis.
#' @param area The area sampled. This is used to express B and P on a per area basis.
#' @param lbl A string that can be used to identify what these calculations are for (perhaps a lake and year label).
#' 
#' @return A list with the following items
#' \describe{
#'   \item{which}{A label for what these calculations were for.}
#'   \item{df}{The original data.frame with mean biomass (\code{mB}), mean weight (\code{mwt}), instantaneous growth rate (\code{G}), and production (\code{P}) columns appended.}
#'   \item{B}{Total biomass for the sampled area (not on a per area basis).}
#'   \item{P}{Total production for the sampled area (not on a per area basis).}
#'   \item{BperA}{Total biomass on a per area basis.}
#'   \item{BperA}{Total biomass on a per area basis.}
#'   \item{PperA}{Total prodouction on a per area per year basis.}
#'   \item{Area}{Area given by the user in \code{area}.}
#' }
#' 
#' @examples 
#' 
calcPB <- function(df,age.c=1,num.c=2,twt.c=3,area=1,adjAgeGaps=TRUE,lbl=NULL) {
  # Some checks
  if (!(is.data.frame(df) | is.matrix(df)))
    stop("'df' must be data.frame or matrix.",call.=FALSE)
  if (is.character(age.c)) age.c <- which(names(df)==age.c)
  if (is.character(num.c)) num.c <- which(names(df)==num.c)
  if (is.character(twt.c)) twt.c <- which(names(df)==twt.c)
  if (length(age.c)<1 | length(num.c)<1 | length(twt.c)<1)
    stop("One of 'age.c', 'num.c', or 'twt.c' is not in 'df'.")
  if (age.c<1 | num.c<1 | twt.c<1)
    stop("One of 'age.c', 'num.c', or 'twt.c' is not a proper column number.")
  if (area<0) stop("'area' must be positive.")
  # Rename columns to match expectations, save old column names
  onames <- names(df)[c(age.c,num.c,twt.c)]
  names(df)[c(age.c,num.c,twt.c)] <- c("age","num","twt")
  # Adjust for gaps in ages
  if (adjAgeGaps) div <- c(NA,diff(df$age))
  else div <- c(NA,rep(1,length(df$age)))
  # Convert total weight to biomass (twt per area)
  df$B <- df$twt/area
  # Add mean biomass (B-bar; div is used to adjust for ages that
  #   represent multiple ages)
  df$mB <- zoo::rollmean(df$B,k=2,na.pad=TRUE,align="right")/div
  # Add mean weight (w-bar)
  df$mwt <- df$twt/df$num
  # Compute instantaneous growth rate (G)
  df$G <- c(NA,diff(log(df$mwt)))
  # Compute age-specific production
  df$P <- df$mB*df$G
  # Put original column names back on data.frame
  names(df)[which(names(df) %in% c("age","num","twt"))] <- onames
  # Create return list
  res <- list(which=lbl,df=as.data.frame(df),
              B=sum(df$twt,na.rm=TRUE)/area,P=sum(df$P,na.rm=TRUE),area=area)
  class(res) <- c("PB","data.frame")
  res
}


plot.PB <- function(PB) {
  ## for the time being, until plotrix is updated
  source("c:/aaaWork/Programs/zWorkingOnR/plotH2.R")
  ## isolate the main data.frame in PB
  df <- PB$df
  par(mfrow=c(1,3))
  par(mar=c(4,4,2,4),mgp=c(1.9,0.5,0),tcl=-0.2,las=1)
  plotH(pnum~age,data=df,xlab="Age (yrs)",
        ylab="Population Numbers",ylim=c(0,1.04*max(df$pnum)),
        yaxs="i")
  par(new=TRUE)
  plot(mwt~age,data=df,type="l",lwd=2,axes=FALSE,xlab=NA,ylab=NA)
  axis(side=4)
  mtext("Mean Weight (kg)",side=4,line=1.9,las=0,cex=par()$cex*par()$cex.axis)
  mtext(paste(PB$which,"-- n =",sum(df$snum),"    PE =",sum(df$pnum)),
        side=3,line=0.5,cex=par()$cex*par()$cex.axis)
  plotH(mB~age,data=df[!is.na(df$mB),],xlab="Age (yrs)",
        ylab="Mean Biomass (kg/ha)",ylim=c(0,1.04*max(df$mB,na.rm=TRUE)),
        yaxs="i")
  par(new=TRUE)
  plot(G~age,data=df[!is.na(df$G),],type="l",lwd=2,axes=FALSE,xlab=NA,ylab=NA)
  axis(side=4)
  mtext("Instantaneous Growth",side=4,line=1.9,las=0,cex=par()$cex*par()$cex.axis)
  mtext(paste("Total Biomass =",formatC(PB$B,format="f",digits=1),"kg/ha"),
        side=3,line=0.5,cex=par()$cex*par()$cex.axis)
  
  par(mar=c(4,6,2,2),mgp=c(1.9,0.5,0),tcl=-0.2,las=1)
  plotH(P~age,data=df,xlab="Age (yrs)",ylab=NA)
  abline(h=0)
  mtext("Production (kg/ha/yr)",side=2,line=2.9,las=0,cex=par()$cex*par()$cex.axis)
  mtext(paste("Total Production =",formatC(PB$P,format="f",digits=2),"kg/ha/yr"),
        side=3,line=0.5,cex=par()$cex*par()$cex.axis)
}