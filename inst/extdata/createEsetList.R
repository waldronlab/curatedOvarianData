#!/usr/bin/env Rscript

##Example invocations of this script:
##
##R CMD BATCH --vanilla "--args patientselection.config patientselection-esets.RData createEsetList.log" createEsetList.R
##
##./createEsetList.R patientselection.config patientselection-esets.RData createEsetList.log

inputArgs <- commandArgs(TRUE)

library(genefilter)
library(survival)
library(logging)


if (length(inputArgs) >=3) {
    kConfigFile <- inputArgs[1]
    source(kConfigFile)
    kOutputFile <- inputArgs[2]
    kLogFile <- inputArgs[3]
} else {
    # assume data was sourced and config options were set
    kLogFile <-"createEsetList.log"
}

## -----------------------------------------------------------------------------
## Logging setup
## -----------------------------------------------------------------------------
basicConfig()
addHandler(writeToFile, logger="", file=kLogFile)

loginfo("Inside script createEsetList.R - inputArgs =")
loginfo(inputArgs)

if (!exists("package.name")) package.name <- "curatedOvarianData"

library(package.name, character.only=TRUE)

loginfo(paste("Loading", package.name, sessionInfo()$otherPkgs[[package.name]]$Version))

## -----------------------------------------------------------------------------
## needed functions
## -----------------------------------------------------------------------------
filterQuantile <- function(object, q){
    if (!identical(q >=0 && q < 1, TRUE))
        stop("require 0 <= q < 1")
    if (!identical(class(object) == "ExpressionSet", TRUE))
        stop("object must be an ExpressionSet")
    gene.sd <- esApply(object,1,sd, na.rm=TRUE)
    gene.quantile <- quantile(gene.sd, probs=q)
    actual.makescutoff <- sum(gene.sd < gene.quantile) / length(gene.sd)
    ##make sure the correct number of genes are getting filtered:
    if (abs(q - actual.makescutoff) > 0.01){
        stop("Not scaling this object, likely pre-scaled.")
    }else{
        object <- object[gene.sd > gene.quantile, ]
    }
    return(object)
}
##recursive intersect function
intersectMany <- function(lst){
    ## Find the intersection of multiple vectors stored as elements of a
    ## list, through a tail-recursive function.
    if (length(lst)==2){
        return(intersect(lst[[1]],lst[[2]]))
    }else{
        return(intersectMany(c(list(intersect(lst[[1]],lst[[2]])),lst[-1:-2])))
    }
}

##Split out non-specific probe sets
expandProbesets <- function (eset, sep = "///"){
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    eset <- eset[order(sapply(x, length)), ]
    x <- lapply(featureNames(eset), function(x) strsplit(x, sep)[[1]])
    idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
    xx <- !duplicated(unlist(x))
    idx <- idx[xx]
    x <- unlist(x)[xx]
    eset <- eset[idx, ]
    featureNames(eset) <- x
    eset
}

## -----------------------------------------------------------------------------
##load the esets
## -----------------------------------------------------------------------------
data(list=data(package=package.name)[[3]][,3])

strEsets <- ls(pattern="^.*_eset$")
esets <- list()

## -----------------------------------------------------------------------------
##Explicit removal of datasets:
## -----------------------------------------------------------------------------
if(exists("remove.datasets") && any(strEsets %in% remove.datasets)){
    remove.datasets <- remove.datasets[remove.datasets %in% strEsets]
    loginfo(paste("Removing ", paste(remove.datasets, collapse=", "), " (remove.datasets)"))
    strEsets <- strEsets[!strEsets %in% remove.datasets]
}

## -----------------------------------------------------------------------------
##Explicit removal of samples from specified datasets:
## -----------------------------------------------------------------------------
delim <- ":"   ##This is the delimiter used to specify dataset:sample,
               ##e.g. TCGA_eset:TCGA.24.1927
if(exists("remove.samples")){
    ##split into those with the delimiter and those without:
    remove.samples.delim <- grep(delim, remove.samples, fixed=TRUE, value=TRUE)
    ##over-write remove.samples for later, keeping this variable in
    ##its historical use:
    remove.samples.orig <- remove.samples
    remove.samples.nodataset <- grep(delim, remove.samples, fixed=TRUE, invert=TRUE, value=TRUE)
    if(length(remove.samples.delim) > 0){
        datasets <- gsub(paste(delim, ".+", sep=""), "", remove.samples.delim)
        samples <- gsub(paste(".+", delim, sep=""), "", remove.samples.delim)
        remove.samples.delim <- lapply(unique(datasets), function(ds){
            samples[datasets %in% ds]
        })
        names(remove.samples.delim) <- unique(datasets)
    }
}

loginfo("Clean up the esets.")
for (strEset in strEsets){
    eset <- get(strEset)
    ##Deal with genes which had a single probe mapping to multiple genes
    if (exists("probes.not.mapped.uniquely")){
        if(identical(probes.not.mapped.uniquely, "drop")){
            ##Drop rows without unique gene name
            eset <- eset[!grepl("///",featureNames(eset),fixed=TRUE),]
        }else if (identical(probes.not.mapped.uniquely, "split")){
            ##Split out rows without unique gene name
            eset <- expandProbesets(eset)
        }
    }
    ## Run ComBat
    if (exists("combat") && combat) {
                                        # workaround bug #12
        if (identical(pubMedIds(eset), "21720365")) {
            eset$batch[match("TCGA.23.1023", sampleNames(eset))] <- 12
        }
        if (sum(is.na(eset$batch)) == 0) {
            tmp <- try(sva::ComBat(exprs(eset),
                                   mod=model.matrix(~rep(1, ncol(eset))), batch=eset$batch), silent=TRUE)
            if (class(tmp)=="matrix") {
                loginfo(paste("Making ComBat correction to", strEset))
                exprs(eset) <- tmp
            }
        }
    }
    ##filter genes with standard deviation below the required quantile
    if(exists("quantile.cutoff") && quantile.cutoff > 0 && quantile.cutoff < 1){
        eset <- filterQuantile(eset, q=quantile.cutoff)
    }
    ##rescale to z-scores
    if(exists("rescale") && rescale){
        exprs(eset) <- t(scale(t(exprs(eset))))
    }
    ##samples to be removed
    remove <- rep(FALSE, ncol(eset))
    ##remove samples without required metadata
    if(exists("meta.required") && length(meta.required) > 0){
        for (varname in meta.required){
            if (varname %in% colnames(pData(eset))){
                remove[ is.na(eset[[varname]]) ] <- TRUE
            }
        }
    }
    ##remove samples not matching regexes
    all.rules <- ls(pattern="rule\\.[0-9]+")
    for (one.rule in all.rules){
        this.remove <- !grepl(get(one.rule)[2], eset[[ get(one.rule)[1] ]])
        if(!strict.checking)
            this.remove[ is.na(eset[[ get(one.rule)[1] ]]) ] <- FALSE
        remove[this.remove] <- TRUE
    }
    ##remove samples pre-specified for removal, that have a dataset specified:
    if(exists("remove.samples.delim")){
        if (strEset %in% names(remove.samples.delim)){
            remove[sampleNames(eset) %in% remove.samples.delim[[strEset]]] <- TRUE
        }
    }
    ##remove samples pre-specified for removal, that did *not* have a dataset specified:
    if(exists("remove.samples.nodataset"))
        remove[sampleNames(eset) %in% remove.samples.nodataset] <- TRUE
    ##do the actual removal
    eset <- eset[, !remove]
    if (exists("considered.datasets") && !(strEset %in% considered.datasets))
    {
        loginfo(paste("excluding",strEset,
                      "(considered.datasets)"))
        next
    }
    ##include study if it has enough samples and events:
    if (exists("min.number.of.events") && !is.na(min.number.of.events)
        && exists("min.sample.size") && !is.na(min.sample.size)
        && min.number.of.events > 0
        && sum(eset$vital_status == "deceased") < min.number.of.events
        || ncol(eset) < min.sample.size)
    {
        loginfo(paste("excluding",strEset,
                      "(min.number.of.events or min.sample.size)"))
        next
    }
    if (exists("min.number.of.genes") && nrow(eset) < min.number.of.genes) {
        loginfo(paste("excluding",strEset,"(min.number.of.genes)"))
        next
    }
    if (exists("remove.retracted") && remove.retracted && length(grep("retracted", experimentData(eset)@other$warnings$warnings)) > 0){
        loginfo(paste("excluding",strEset,"(remove.retracted)"))
        next
    }
    if (exists("remove.subsets") && remove.subsets && length(grep("subset", experimentData(eset)@other$warnings$warnings)) > 0){
        loginfo(paste("excluding",strEset,"(remove.subsets)"))
        next
    }
    loginfo(paste("including",strEset))
    ##    featureNames(eset) <- make.names(featureNames(eset))  ##should not do this, it is irreversible.
    esets[[strEset]] <- eset
    rm(eset)
}

##optionally take the intersection of genes common to all platforms:
if(exists("keep.common.only") && keep.common.only){
    features.per.dataset <- lapply(esets, featureNames)
    intersect.genes <- intersectMany(features.per.dataset)
    esets <- lapply(esets, function(eset){
        eset <- eset[intersect.genes, ]
        return(eset)
    })
}

ids.with.missing.data <- which(sapply(esets, function(X)
                                      sum(!complete.cases(exprs(X))) > 0))
loginfo(paste("Ids with missing data:", paste(names(ids.with.missing.data),
                                              collapse=", ")))

if (length(ids.with.missing.data) > 0 && exists("impute.missing") && impute.missing) {
    for (i in ids.with.missing.data) {
        require(impute)
        exprs(esets[[i]]) = impute.knn(exprs(esets[[i]]))$data
    }
}


if (exists("add.surv.y") && is.function(add.surv.y)) {
    for (i in 1:length(esets)) {
        esets[[i]]$y = add.surv.y(esets[[i]])
    }
}


if (exists("kOutputFile"))
    save(esets,file=kOutputFile)
