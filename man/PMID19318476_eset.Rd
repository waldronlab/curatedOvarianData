\name{PMID19318476_eset}
\alias{PMID19318476_eset}
\docType{data}
\title{Microarray analysis of early stage serous ovarian cancers shows profiles predictive of favorable outcome.}
\description{Although few women with advanced serous ovarian cancer are cured, detection of the disease at an early stage is associated with a much higher likelihood of survival. We previously used gene expression array analysis to distinguish subsets of advanced cancers based on disease outcome. In the present study, we report on gene expression of early-stage cancers and validate our prognostic model for advanced-stage cancers.}
\usage{data(PMID19318476_eset)}
\format{

Authors: Berchuck A, Iversen ES, Luo J, Clarke JP, Horne H, Levine DA, Boyd J, Alonso MA, Secord AA, Bernardini MQ, Barnett JC, Boren T, Murphy SK, Dressman HK, Marks JR, Lancaster JM.

Lab: Department of Obstetrics and Gynecology, Division of Gynecologic Oncology, Institute for Genome Sciences and Policy, Duke University Medical Center, Durham, North Carolina 27710, USA.

Contact: berch001@mc.duke.edu

Title: Microarray analysis of early stage serous ovarian cancers shows profiles predictive of favorable outcome.

Url: http://pubget.com/search?q=19318476

PubMedID: 19318476
}
\details{
Platform used: HG-U133A.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 42 samples.

GEO_platform_accession: PMID19318476
             
platform_summary: GPL96
             
biomart_ID: HG-U133A affy_hg_u133a

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

records = 42.00  n.max = 42.00   n.start = 42.00   events = 22.00  median = 2.79    0.95LCL = 2.30 0.95UCL = NA 

--------------------------------------------------------------------------------- 


Available sample meta-data: 

---------------------------------------------------------------------------------

alt_sample_name:
Length = 42.
Class = character.
Mode = character.

summarygrade: 
high = 24. low = 17. NA's = 1.

summarystage: 
early = 2. late = 39. NA's = 1.

tumorstage:
1 = 1.  2  = 1.  3 = 29. 4 = 10. NA's = 1.

substage: 
a = 1. b = 1. c = 29. NA's = 11.

G: 
1 = 2.  2 = 15.  3 = 24. NA's = 1.
 
age_at_initial_pathologic_diagnosis: 
Min. = 33.00  1st Qu. = 55.00 Median =  62.00 Mean = 61.46  3rd Qu. = 70.00   Max. = 81.00 NA's: 1.00 

recurrence_status: 
norecurrence = 6. recurrence = 36.

days_to_death: 
Min. = 30.0  1st Qu. = 367.5  Median = 825.0  Mean = 1105.0  3rd Qu. = 1050.0  Max. = 3420.0 

vital_status: 
deceased = 22.  living = 20.

debulking: 
optimal = 20. suboptimal = 21. NA's = 1.

batch: 

2004-03-09 2004-03-16 2004-04-20 2004-05-18 2004-05-21 2004-05-27 2004-06-22 

        14          3          4          8          6          5          1 
        
2004-06-23 

         1 
}
\examples{
data(PMID19318476_eset)
## maybe str(PMID19318476_eset) ; plot(PMID19318476_eset) ...

if(require(affy)){
  summary(PMID19318476_eset$recurrence_status)
  time <- PMID19318476_eset$days_to_death / 365
  cens <- ifelse(PMID19318476_eset$vital_status=="deceased",1,0)
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
