\name{GSE19829.GPL570_eset}
\alias{GSE19829.GPL570_eset}
\docType{data}
\title{Gene Expression Profile of BRCAness That Correlates With Responsiveness to Chemotherapy and With Outcome in Patients With Epithelial Ovarian Cancer.}
\description{To define a gene expression profile of BRCAness that correlates with chemotherapy response and outcome in epithelial ovarian cancer (EOC).}
\usage{data(GSE19829.GPL570_eset)}
\format{

Authors: Konstantinopoulos PA, Spentzos D, Karlan BY, Taniguchi T et al.

Lab: Beth Israel Deaconess Medical Center, Harvard Medical School, Boston, MA, USA.

Contact: scannist@bidmc.harvard.edu

Title: Gene Expression Profile of BRCAness That Correlates With Responsiveness to Chemotherapy and With Outcome in Patients With Epithelial Ovarian Cancer.

Url: http://pubget.com/search?q=20547991

PubMedID: 20547991


}
\details{
Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.

frma normalization used from the frma bioconductor package.

assayData: 18708 features, 28 samples.

GEO_platform_accession: GSE19829-GPL570.
platform_summary: GPL570
biomart_ID: HG-U133_Plus_2 affy_hg_u133_plus_2

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

records = 28.00 n.max = 28.00  n.start = 28.00  events = 17.00  median =  3.95    0.95LCL = 2.71 0.95UCL =  NA 

------------------------------------------------------------------------------ 


Available sample meta-data: 

------------------------------------------------------------------------------

alt_sample_name: 
Length = 28.
Class = character.
Mode = character. 

sample_type: 
tumor = 28.
   
primarysite: 
ov = 28.

days_to_death: 
Min. = 150 1st Qu.= 540 Median = 1050 Mean = 1291 3rd. Qu = 1688 Max. = 3450 

vital_status: 
deceased = 17. living = 11.
      
debulking: 
unknown = 28.
    

batch: 

2009-08-14 

        28 

}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE19829
}
\examples{
data(GSE19829.GPL570_eset)
## maybe str(GSE19829-GPL570_eset) ; plot(GSE19829-GPL570_eset) ...

if(require(affy)){
  summary(GSE19829.GPL570_eset$vital_status)
  time <- GSE19829.GPL570_eset$days_to_death / 365
  cens <- ifelse(GSE19829.GPL570_eset$vital_status=="deceased",1,0)
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
