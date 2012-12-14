\name{GSE14764_eset}
\alias{GSE14764_eset}
\docType{data}
\title{A prognostic gene expression index in ovarian cancer - validation across different independent data sets.}
\description{Ovarian carcinoma has the highest mortality rate among gynaecological malignancies. In this project, we investigated the hypothesis that molecular markers are able to predict outcome of ovarian cancer independently of classical clinical predictors, and that these molecular markers can be validated using independent data sets. We applied a semi-supervised method for prediction of patient survival. Microarrays from a cohort of 80 ovarian carcinomas (TOC cohort) were used for the development of a predictive model, which was then evaluated in an entirely independent cohort of 118 carcinomas (Duke cohort). A 300-gene ovarian prognostic index (OPI) was generated and validated in a leave one out approach in the TOC cohort (Kaplan-Meier analysis, p = 0.0087). In a second validation step, the prognostic power of the OPI was confirmed in an independent data set (Duke cohort, p = 0.0063). In multivariate analysis, the OPI was independent of the post-operative residual tumour, the main clinico-pathological prognostic parameter with an adjusted hazard ratio of 6.4 (TOC cohort, CI 1.8-23.5, p = 0.0049) and 1.9 (Duke cohort, CI 1.2-3.0, p = 0.0068). We constructed a combined score of molecular data (OPI) and clinical parameters (residual tumour), which was able to define patient groups with highly significant differences in survival. The integrated analysis of gene expression data as well as residual tumour can be used for optimized assessment of the prognosis of platinum taxol treated ovarian cancer. As traditional treatment options are limited, this analysis may be able to optimize clinical management and to identify those patients who would be candidates for new therapeutic strategies.}
\usage{data(GSE14764_eset)}
\format{

Authors: Denkert C, Budczies J, Darb-Esfahani S, Gyorffy B et al.

Lab: Institute of Pathology, Charite University Hospital, Berlin, Germany.

Contact: carsten.denkert@charite.de

Title: A prognostic gene expression index in ovarian cancer: Validation across different independent data sets.

Url: http://pubget.com/search?q=19294737

PubMedID: 19294737
}
\details{
Platform used: [HG-U133A] Affymetrix Human Genome U133A Array.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 80 samples.

GEO_platform_accession: GSE14764.
platform_summary: GPL96.
biomart_ID: HG-U133A affy_hg_u133a.

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

records = 80.00   n.max = 80.00   n.start = 80.00   events = 21.00  median = 4.52 0.95LCL = 4.19    0.95UCL = NA 

---------------------------------------------------------------------------------

Available sample meta-data: 


---------------------------------------------------------------------------------
alt_sample_name: 
Min. = 1.00   1st Qu. = 20.75  Median = 40.50   Mean = 40.50   3rd Qu. = 60.25   Max. = 80.00 

sample_type: 
tumor = 80.

subtype: 
Length = 80.
Class = character.
Mode = character.

primarysite: 
ov = 80.

summarygrade: 
high = 54.  low = 26.

summarystage: 
early = 9. late = 71.

tumorstage:
1 = 8. 2 = 1. 3 = 69. 4 = 2.

substage: 
a = 4. b  = 6. c  = 32. NA's = 38.

G: 
1 = 3. 2 = 23. 3 = 54.

recurrence_status: 
norecurrence = 50.  recurrence = 26. NA's = 4.

days_to_death: 
Min. = 210. 1st Qu. = 660. Median = 1050. Mean = 1011. 3rd Qu. = 1328. Max = 2190.

vital_status:
deceased = 21. living = 59.

debulking: 
unknown = 80.

batch: 

2004-09-29 2004-09-30 2004-10-01 2005-01-21 2005-01-25 2005-01-26 2005-01-28 

         1          2          6          4          7          8         10 
         
2005-03-02 2006-07-26 2006-07-27 2006-07-28 2006-08-11 2006-08-18 2006-08-19 

         6          4          6          4         10          3          4 
         
2006-08-21 

         5 
        
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE14764
}
\examples{
data(GSE14764_eset)
## maybe str(GSE14764_eset) ; plot(GSE14764_eset) ...
if(require(affy)){
  summary(GSE14764_eset$recurrence_status)
  time <- GSE14764_eset$days_to_death / 365
  cens <- ifelse(GSE14764_eset$vital_status=="deceased",1,0)
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
