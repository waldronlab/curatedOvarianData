\name{GSE30009_eset}
\alias{GSE30009_eset}
\docType{data}
\title{Multidrug Resistance-Linked Gene Signature Predicts Overall Survival of Patients with Primary Ovarian Serous Carcinoma}
\description{qRT-PCR dataset.  This study assesses the ability of multidrug resistance (MDR)-associated gene expression patterns to predict survival in patients with newly diagnosed carcinoma of the ovary. The scope of this research differs substantially from that of previous reports, as a very large set of genes was evaluated whose expression has been shown to affect response to chemotherapy.}

\usage{data(GSE30009_eset)}
\format{

Note this is a qRT-PCR dataset, not microarray.
  
Authors:  	Gillet J, Calcagno AM, Varma S, Davidson B, Elstrand MB, Ganapathi R, Kamat A, Sood A, Ambudkar SV, Seiden M, Rueda B, Gottesman MM

Contact: Jean-Pierre Gillet <gilletjp@mail.nih.gov>

Title: Multidrug resistance-linked gene signature predicts overall survival of patients with primary ovarian serous carcinoma.

Url: http://pubget.com/search?q=22492981

PubMedID: 22492981
}
\details{
  \preformatted{
Platform used: TaqMan qRT-PCR Homo sapiens Low-Density Array 380

Processed data here are those provided by the original study authors.

assayData: 380 features, 103 samples.

GEO_platform_accession: GPL13728


---------------------------------------------------------------------------------

Available sample meta-data: 

---------------------------------------------------------------------------------

sample_type
tumor 
  103 

subtype
clearcell       ser 
        1       102 

summarygrade
high  low 
  92    9 

summarystage
late 
 103 

T   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  3.000   3.000   3.000   3.204   3.000   4.000 

substage
 b  c 
 2 60 

G   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  1.000   3.000   3.000   2.871   3.000   3.000       2 

age_at_initial_pathologic_diagnosis   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  30.00   56.00   61.00   62.45   71.50   87.00 

days_to_death   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     24     598    1053    1156    1568    4748 

vital_status
deceased   living 
      57       46 

debulking
   optimal suboptimal 
   81         22
   }
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE30009
}
\examples{
data(GSE30009_eset)
meta.data <- pData(GSE30009_eset)
for (i in 1:ncol(meta.data)){
    if (!all(is.na(meta.data[[i]])) & length(unique(meta.data[[i]])) < length(meta.data[[i]])){
        cat(colnames(meta.data)[i])
        if(class(meta.data[[i]]) == "numeric" || class(meta.data[[i]]) == "integer"){
            print(summary(meta.data[[i]]))
        }else if(class(meta.data[[i]]) == "character" || class(meta.data[[i]]) == "factor" || class(meta.data[[i]]) == "logical"){
            print(table(meta.data[[i]]))
        }
        cat("\n")
    }
}
if(require(affy)){
  summary(GSE30009_eset$vital_status)
  time <- GSE30009_eset$days_to_death / 365
  cens <- ifelse(GSE30009_eset$vital_status=="deceased",1,0)
  library(survival)
  fit <- survfit(Surv(time,cens)~-1)
  fit          #for summary of survival
  summary(fit)
  plot(fit,xlab="Time (years)",ylab="Survivor function")
  inverse.fit <- survfit(Surv(time,1-cens) ~ -1)
  inverse.fit  #for summary of follow-up time
}
}
\keyword{datasets}
