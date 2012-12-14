\name{PMID17290060_eset}
\alias{PMID17290060_eset}
\docType{data}
\title{An integrated genomic-based approach to individualized treatment of patients with advanced-stage ovarian cancer.}
\description{The purpose of this study was to develop an integrated genomic-based approach to personalized treatment of patients with advanced-stage ovarian cancer. We have used gene expression profiles to identify patients likely to be resistant to primary platinum-based chemotherapy and also to identify alternate targeted therapeutic options for patients with de novo platinum-resistant disease.}
\usage{data(PMID17290060_eset)}
\format{

Authors: Dressman HK, Berchuck A, Chan G, Zhai J, Bild A, Sayer R, Cragun J, Clarke J, Whitaker RS, Li L, Gray J, Marks J, Ginsburg GS, Potti A, West M, Nevins JR, Lancaster JM.

Lab: Division of Gynecologic Surgical Oncology, H. Lee Moffitt Cancer Center and Research Institute, Tampa, FL 33612, USA.

Contact: lancasjm@moffitt.usf.edu

Title: An integrated genomic-based approach to individualized treatment of patients with advanced-stage ovarian cancer.

Url: http://pubget.com/search?q=17290060

PubMedID: 17290060

}
\details{
Platform used: HG-U133A.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 117 samples.

GEO_platform_accession: PMID17290060
platform_summary: GPL96
biomart_ID: HG-U133A affy_hg_u133a

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

records = 117.00  n.max = 117.00  n.start = 117.00  events = 67.00 median = 5.26    0.95LCL = 2.79 0.95UCL = 7.48

--------------------------------------------------------------------------------


Available sample meta-data: 

--------------------------------------------------------------------------------

alt_sample_name: 
Length = 117.
Class = character.
Mode = character.

primarysite: 
ov = 117.

summarygrade: 
high = 57. low = 57. NA's = 3.

summarystage: 
early = 1. late = 115. NA's = 1.

tumorstage:
2 = 1. 3 = 98.  4 = 17. NA's = 1.

G: 
1 = 4. 2  = 53. 3 = 56. 4 = 1. NA's = 3.

days_to_death: 
Min. = 30. 1st Qu. = 510. Median = 1020. Mean =  1496. 3rd Qu. = 2220. Max. = 5550.

vital_status: 
deceased = 67.  living = 50.

primary_therapy_outcome_success: 
completeresponse = 85. progressivedisease = 32.

debulking: 
optimal = 63. suboptimal = 54.
}
\examples{
data(PMID17290060_eset)
## maybe str(PMID17290060_eset) ; plot(PMID17290060_eset) ...

##ISSUE: only "primary_therapy_outcome_success"

if(require(affy)){
  summary(PMID17290060_eset$primary_therapy_outcome_success)
  time <- PMID17290060_eset$days_to_death / 365
  cens <- ifelse(PMID17290060_eset$vital_status=="deceased",1,0)
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
