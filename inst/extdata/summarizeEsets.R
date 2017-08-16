esetToSurv <- function(eset, output, months=TRUE){
  if (sum(eset$vital_status =="deceased", na.rm=TRUE) == 0) return(NULL)
  eset <- eset[ ,!is.na(eset$vital_status) & !is.na(eset$days_to_death)]
  time <- eset$days_to_death
  if (months)
      time <- time / 30.4
  cens <- rep(NA,length(time))
  if(output == "event"){
      cens[eset$vital_status == "living"] <- 0
      cens[eset$vital_status == "deceased"] <- 1
  }else if(output == "reverseevent"){
      cens[eset$vital_status == "living"] <- 1
      cens[eset$vital_status == "deceased"] <- 0
  }
  Surv(time,cens)
}
medianSurvival <- function(eset){
    ##make surv objects
    y <- esetToSurv(eset, "event")
    if (is.null(y)) {
        binary.table <- paste(table(eset$os_binary), collapse="/")
        if(length(binary.table) == 0 || binary.table == "")
            binary.table <- NA
        output <- c(rep(NA,3), binary.table )
    } else {
        y.rev <- esetToSurv(eset, "reverseevent")
        ##kaplan-meier fits
        km.fit <- survfit(y ~ 1)
        reverse.km.fit <- survfit(y.rev ~ 1)
        ##get median survival and follow-up times
        surv.median.position <- min( which(km.fit$surv < 0.5) )
        surv.median.time <- km.fit$time[ surv.median.position ]
        followup.median.position <- min( which(reverse.km.fit$surv < 0.5) )
        followup.median.time <- reverse.km.fit$time[ followup.median.position ]
        ## % censoring
        percent.censoring <- round(sum(y[,2] == 0) / nrow(y) * 100)
        output <- round( c(surv.median.time, followup.median.time, percent.censoring, NA) )
    }
    names(output) <- c("median.survival", "median follow-up", "percent censoring", "binarized OS (long/short)")
    return(output)
}
getEsetData <- function(eset, 
    possible.stages = c("early", "late"),
    possible.histological_types = c("ser", "clearcell", "endo", "mucinous", "other", "unknown")
){
    makeSummary <- function(field, possible.values){
        my.table <- sapply(possible.values, function(x) sum(eset[[field]] == x, na.rm=TRUE))
        my.table <- c(my.table, ncol(eset) - sum(my.table))
        my.table <- paste(my.table, collapse="/")
        return(my.table)
    }
    survdata <- medianSurvival(eset)
    output <- c(experimentData(eset)@lab,  #name/year
                experimentData(eset)@pubMedIds,
                ncol(eset),  # number of samples
                makeSummary("summarystage", possible.stages),
                makeSummary("histological_type", possible.histological_types),
                sum(eset$histological_type == possible.histological_types[1] &
                eset$summarystage == possible.stages[2], na.rm=TRUE), ##number of serous samples
                ifelse(length(experimentData(eset)@other$platform_shorttitle), experimentData(eset)@other$platform_shorttitle,""),  #platform
                survdata)
    names(output) <- c("Study",
                       "PMID",
                       "N samples",
                       "stage",
                       "histology",
                       "known late-stage serous",
                       "Platform",
                       names(survdata))
    return(output)
}
