\name{GSE19829.GPL8300_eset}
\alias{GSE19829.GPL8300_eset}
\docType{data}
\title{Gene Expression Profile of BRCAness That Correlates With Responsiveness to Chemotherapy and With Outcome in Patients With Epithelial Ovarian Cancer.}
\description{To define a gene expression profile of BRCAness that correlates with chemotherapy response and outcome in epithelial ovarian cancer (EOC).}
\usage{data(GSE19829.GPL8300_eset)}
\format{

Authors: Konstantinopoulos PA, Spentzos D, Karlan BY, Taniguchi T et al.

Lab: Beth Israel Deaconess Medical Center, Harvard Medical School, Boston, MA, USA.

Contact: scannist@bidmc.harvard.edu

Title: Gene Expression Profile of BRCAness That Correlates With Responsiveness to Chemotherapy and With Outcome in Patients With Epithelial Ovarian Cancer.

Url: http://pubget.com/search?q=20547991

PubMedID: 20547991

}
\details{
Platform used: [HG_U95Av2] Affymetrix Human Genome U95 Version 2 Array.

rma normalization used from the rma bioconductor package.

assayData: 8913 features, 42 samples.

GEO_platform_accession: GSE19829-GPL8300
platform_summary: GPL8300
biomart_ID: HG_U95Av2 affy_hg_u95av2

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

records = 42.00 n.max = 42.00 n.start = 42.00 events = 23.00 median = 3.78    0.95LCL = 2.79 0.95UCL =  NA 

---------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------

alt_sample_name: 
Length = 42.
Class = character.
Mode = character.

sample_type: 
tumor = 42.

primarysite: 
ov = 42.

days_to_death: 
Min. = 30.0 1st Qu. = 727.5 Median = 1155.0 Mean = 1089.0 3rd Qu. = 1485.0  Max. = 2040.0

vital_status: 
deceased = 23. living = 19.

debulking: 
unknown = 42.

batch: 

2001-09-14 2001-12-14 2002-08-20 2003-09-09 2003-09-18 

         7          4         14         13          4 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE19829
}

\examples{
data(GSE19829.GPL8300_eset)
## maybe str(GSE19829-GPL8300_eset) ; plot(GSE19829-GPL8300_eset) ...
if(require(affy)){
  summary(GSE19829.GPL8300_eset$vital_status)
  time <- GSE19829.GPL8300_eset$days_to_death / 365
  cens <- ifelse(GSE19829.GPL8300_eset$vital_status=="deceased",1,0)
  library(survival)
  fit <- survfit(Surv(time,cens)~-1)
  ##fit          #for summary of survival
  summary(fit)
  plot(fit,xlab="Time (years)",ylab="Survivor function")
  inverse.fit <- survfit(Surv(time,1-cens) ~ -1)
  inverse.fit  #for summary of follow-up time
 }
}
\keyword{datasets}
