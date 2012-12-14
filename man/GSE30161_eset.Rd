\name{GSE30161_eset}
\alias{GSE30161_eset}
\docType{data}
\title{Multi-gene expression predictors of single drug responses to adjuvant chemotherapy in ovarian carcinoma: predicting platinum resistance}
\description{We adapted and applied the previously-published COXEN algorithm to develop molecular predictors for therapeutic responses of patients' tumors based on expression signatures derived from the NCI-60 in vitro drug activities and genomic expression data. Genome-wide candidate biomarkers were first triaged by examining expression patterns of frozen and formalin-fixed paraffin embedded (FFPE) tissue samples. We then identify initial drug sensitivity biomarkers for carboplatin and paclitaxel, respectively. These biomarkers were further narrowed by examining concordant expression patterns between cell lines and a historical set of ovarian cancer patients. Multivariate predictors were obtained from the NCI-60 cell lines and refined using historical patient cohorts. To independent validate these molecular predictors, we performed genome-wide profiling on FFPE samples of 58 ovarian cancer patients obtained prior to adjuvant chemotherapy.}
\usage{data(GSE30161_eset)}
\format{

Authors: Ferriss JS, Kim Y, Duska L, Birrer M, Levine DA, Moskaluk C, Theodorescu D, Lee JK.

Lab: University of Virginia, Department of Public Health Sciences

Contact: Youngchul Kim

Title: Multi-gene expression predictors of single drug responses to adjuvant chemotherapy in ovarian carcinoma: predicting platinum resistance.

PubMedID: 22348014

}
\details{
\preformatted{
Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.

frma normalization used from the frma bioconductor package.

assayData: 18708 features, 63 samples.

GEO_platform_accession: GSE30161

platform_summary: GPL570

biomart_ID: HG-U133_Plus_2 affy_hg_u133_plus_2

assayData: 19093 features, 58 samples

         GEO_platform_accession platform_summary          biomart_ID
GSE30161                 GPL570   HG-U133_Plus_2 affy_hg_u133_plus_2

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  58.00   58.00   58.00   36.00    4.19    2.70    6.17 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       58 character character 

sample_type: 
tumor 
   58 

subtype: 
   Length     Class      Mode 
       58 character character 

summarygrade: 
high  low NA's 
  33   21    4 

summarystage: 
late 
  58 

tumorstage: 
 3  4 
53  5 

substage: 
 a  b  c 
 9 11 38 

G: 
   1    2    3 NA's 
   2   19   33    4 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  38.00   53.50   62.00   62.57   72.00   85.00 

pltx: 
 y 
58 

tax: 
 n  y 
 4 54 

neo: 
 n 
58 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   12.0   255.2   386.0   742.1   768.2  4208.0 

recurrence_status: 
norecurrence   recurrence         NA's 
           6           48            4 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   49.0   585.2  1010.0  1375.0  2131.0  4208.0 

vital_status: 
deceased   living 
      36       22 

debulking: 
   optimal suboptimal       NA's 
        26         30          2 
}

}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE30161
}
\examples{
data(GSE30161_eset)
## maybe str(GSE30161_eset) ; plot(GSE30161_eset) ...
if(require(affy)){
  summary(GSE30161_eset$vital_status)
  time <- GSE30161_eset$days_to_death / 365
  cens <- ifelse(GSE30161_eset$vital_status=="deceased",1,0)
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
