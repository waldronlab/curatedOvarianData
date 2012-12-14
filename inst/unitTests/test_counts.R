test_counts <- function() {
    require(curatedOvarianData)
    require(affy)
    strEsets <- data(package="curatedOvarianData")[[3]][,3]
    esets <- lapply(strEsets, function(x){
        data(list=x)
        get(x)
    })
    df <- data.frame(PMID=sapply(esets, pubMedIds), 
                     ncol=sapply(esets,ncol),
                     nrow=sapply(esets,nrow),
                     nwithbatch=sapply(esets, function(x) sum(!is.na(x$batch))),
                     deceased=sapply(esets,function(X)
                               sum(X$vital_status=="deceased",na.rm=TRUE)),
                     stringsAsFactors = FALSE)
    df = df[order(df$PMID, df$ncol, df$nrow),]
    ## reference file was generated with the following commands on 0.9
    ##  write.csv(df, file="../extdata/curatedOvarianData_counts.csv", quote=TRUE, row.names=FALSE)
    ## dfref <- read.csv("../extdata/curatedOvarianData_counts.csv", as.is=TRUE)
    dfref <- read.csv(system.file("extdata/curatedOvarianData_counts.csv", package = "curatedOvarianData"), as.is=TRUE)
    sapply(1:ncol(df), function(i){
        df[df[, i] != dfref[ ,i], i] <- paste("!!!", df[df[, i] != dfref[ ,i], i], "!!!", sep="")
        print(data.frame(i, df[, i], dfref[, i]))
        checkEquals(as.character(df[,i]), as.character(dfref[,i]))})
}
