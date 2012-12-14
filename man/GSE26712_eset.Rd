\name{GSE26712_eset}
\alias{GSE26712_eset}
\docType{data}
\title{A gene signature predicting for survival in suboptimally debulked patients with ovarian cancer.}
\description{Despite the existence of morphologically indistinguishable disease, patients with advanced ovarian tumors display a broad range of survival end points. We hypothesize that gene expression profiling can identify a prognostic signature accounting for these distinct clinical outcomes. To resolve survival-associated loci, gene expression profiling was completed for an extensive set of 185 (90 optimal/95 suboptimal) primary ovarian tumors using the Affymetrix human U133A microarray. Cox regression analysis identified probe sets associated with survival in optimally and suboptimally debulked tumor sets at a P value of <0.01. Leave-one-out cross-validation was applied to each tumor cohort and confirmed by a permutation test. External validation was conducted by applying the gene signature to a publicly available array database of expression profiles of advanced stage suboptimally debulked tumors. The prognostic signature successfully classified the tumors according to survival for suboptimally (P = 0.0179) but not optimally debulked (P = 0.144) patients. The suboptimal gene signature was validated using the independent set of tumors (odds ratio, 8.75; P = 0.0146). To elucidate signaling events amenable to therapeutic intervention in suboptimally debulked patients, pathway analysis was completed for the top 57 survival-associated probe sets. For suboptimally debulked patients, confirmation of the predictive gene signature supports the existence of a clinically relevant predictor, as well as the possibility of novel therapeutic opportunities. Ultimately, the prognostic classifier defined for suboptimally debulked tumors may aid in the classification and enhancement of patient outcome for this high-risk population.}
\usage{data(GSE26712_eset)}
\format{

Authors: Bonome T, Levine DA, Shih J, Randonovich M et al

Lab: Cell and Cancer Biology Branch, National Cancer Institute, NIH, Rockville, Maryland 20892, USA.

Contact: birrerm@mail.nih.gov

Title: A gene signature predicting for survival in suboptimally debulked patients with ovarian cancer.

Url: http://pubget.com/search?q=18593951

PubMedID: 18593951

 
}
\details{
Platform used: [HG-U133A] Affymetrix Human Genome U133A Array.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 195 samples.

GEO_platform_accession: GSE26712

platform_summary: GPL96

biomart_ID: HG-U133A affy_hg_u133a

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

10 observations deleted due to missingness 

records = 185.00  n.max = 185.00  n.start = 185.00  events = 129.00  median = 3.83    0.95LCL = 3.24   0.95UCL = 4.83 

--------------------------------------------------------------------------------


Available sample meta-data: 

--------------------------------------------------------------------------------

alt_sample_name: 
Length = 195.
Class = character.
Mode = character.
     
sample_type: 
healthy = 10.  tumor = 185.

subtype: 
ser = 185. NA's = 10.

primarysite: 
ov = 195.

summarygrade: 
high = 185. NA's = 10.

summarystage: 
late = 185. NA's = 10.

tumorstage:
4 = 185. NA's = 10.

G: 
4 = 185. NA's = 10.

recurrence_status: 
norecurrence = 42. recurrence = 153.

days_to_death: 
Min. = 21.9   1st Qu. = 660.6 Median = 1164.0 Mean = 1429.0 3rd Qu. = 1880.0  Max. = 4982.0 NA's = 10.0 

vital_status: 
deceased = 129.  living = 56.  NA's = 10.

debulking: 
optimal = 90. suboptimal = 95. unknown = 10.

batch: 

2003-11-04 2003-11-05 2003-11-06 2003-11-07 2003-11-20 2003-11-21 2003-12-16 

        14         16          9          6         10         15         17 
        
2003-12-23 2003-12-24 2004-04-20 2004-04-21 2004-04-27 2004-09-28 2005-07-27 

        12         11         20         17          9         14         15 
        
2006-11-09 

        10 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE26712
}
\examples{
data(GSE26712_eset)
## maybe str(GSE26712_eset) ; plot(GSE26712_eset) ...

if(require(affy)){
  summary(GSE26712_eset$recurrence_status)
  time <- GSE26712_eset$days_to_death / 365
  cens <- ifelse(GSE26712_eset$vital_status=="deceased",1,0)
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
