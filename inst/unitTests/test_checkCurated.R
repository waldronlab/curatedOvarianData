test_checkCurated <- function(){
    library(affy)
    ##the template to use
    ##template.file <- "../extdata/template_ov.csv"
    template.file <- system.file("extdata/template_ov.csv", package = "curatedOvarianData")
    template <- read.csv(template.file, as.is=TRUE)

    ##get rid of the first row of template - this is sample_name
    template <- template[-1, ]

    ##the list of datasets to check
    dataset.names <- data(package="curatedOvarianData")$results[, "Item"]

    ##Create a checklist for three checks in all files
    checklist <- matrix(NA,ncol=2,nrow=length(dataset.names))
    colnames(checklist) <- c("column.names","variable.regex.checks")
    rownames(checklist) <- dataset.names

    ##Make a list of curated metadata:
    curated <- lapply(dataset.names, function(x){
        data(list=x, package="curatedOvarianData")
        pData(get(x))
    })
    names(curated) <- dataset.names


    ##-------------------------------------------------
    ##check that the column names match the template
    ##-------------------------------------------------
    for (i in 1:length(curated)){
        print(paste("Checking column names in",dataset.names[i]))
        ##check column names
        if(identical(all.equal(colnames(curated[[i]]),template$col.name),TRUE)){
            print("column names OK")
            checklist[i,"column.names"] <- TRUE
        }else{
            all.equal(colnames(curated[[i]]),template$col.name)
            print("Column names should be: ")
            print(template$col.name)
            print("Column names are: ")
            print(colnames(curated[[i]]))
            checklist[i,"column.names"] <- FALSE
        }
    }

    ##-------------------------------------------------
    ##construct the regexes from template$allowedvalues
    ##-------------------------------------------------
    regexes <- template$allowedvalues
    regexes <- paste("^",regexes,"$",sep="")
    regexes <- gsub("|","$|^",regexes,fixed=TRUE)
    ##regexes[template$requiredness=="optional"] <- paste(regexes[template$requiredness=="optional"],"|^NA$",sep="")
    regexes[grep("^*$",regexes,fixed=TRUE)] <- "."
    names(regexes) <- template$col.name
    regexes[grep("decimal",regexes)] <- "^[0-9]+\\.?[0-9]*$"

    ##-------------------------------------------------
    ##Check the data entries in each column for regex
    ## matching, uniqueness, and missingness
    ##-------------------------------------------------
    for (i in which(checklist[,"column.names"])){
        print(paste("Checking regexes in",dataset.names[i]))
        ##check column names
        column.OK <- rep(NA,ncol(curated[[i]]))
        names(column.OK) <- colnames(curated[[i]])
        for(j in 1:ncol(curated[[i]])){
            doesmatch <- grep(regexes[j],curated[[i]][,j])
            if(template[j,"requiredness"]=="optional"){
                doesmatch <- c(doesmatch,which(is.na(curated[[i]][,j])))
            }
            doesnotmatch <- 1:nrow(curated[[i]])
            doesnotmatch <- doesnotmatch[!doesnotmatch %in% doesmatch]
            ## if field must be unique, add non-unique values to doesnotmatch
            if(template[j,"uniqueness"]=="unique"){
                counts.table <- table(curated[[i]][,j])
                counts <- counts.table[match(curated[[i]][,j],names(counts.table))]  #this counts the occurences
                counts[is.na(counts)] <- 1  ##Consider NAs to be unique
                nonunique <- which(counts>1)
                doesnotmatch <- c(doesnotmatch,nonunique)
                doesnotmatch <- unique(doesnotmatch)  #don't duplicate fields which fail both tests
            }
            if(length(doesnotmatch)==0){
                column.OK[j] <- TRUE
            }else{
                column.OK[j] <- FALSE
                curated[[i]][doesnotmatch,j] <- paste("!!!",curated[[i]][doesnotmatch,j],"!!!",sep="")
            }
        }
        if(all(column.OK)){
            print("all regex checks OK")
            checklist[i,"variable.regex.checks"] <- TRUE
        }else{
            print(paste("The following columns failed the regex check for ",dataset.names[i],":",sep=""))
            print(names(column.OK)[!column.OK])
            print(curated[[i]][,!column.OK])
            checklist[i,"variable.regex.checks"] <- FALSE
        }
    }

    print(checklist)

    checkTrue( all(checklist[, "column.names"]) )
    checkTrue( all(checklist[, "variable.regex.checks"]) )

}
