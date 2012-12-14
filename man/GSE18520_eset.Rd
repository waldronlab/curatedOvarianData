\name{GSE18520_eset}
\alias{GSE18520_eset}
\docType{data}
\title{A gene signature predictive for outcome in advanced ovarian cancer identifies a survival factor: microfibril-associated glycoprotein 2.}
\description{Advanced stage papillary serous tumors of the ovary are responsible for the majority of ovarian cancer deaths, yet the molecular determinants modulating patient survival are poorly characterized. Here, we identify and validate a prognostic gene expression signature correlating with survival in a series of microdissected serous ovarian tumors. Independent evaluation confirmed the association of a prognostic gene microfibril-associated glycoprotein 2 (MAGP2) with poor prognosis, whereas in vitro mechanistic analyses demonstrated its ability to prolong tumor cell survival and stimulate endothelial cell motility and survival via the alpha(V)beta(3) integrin receptor. Increased MAGP2 expression correlated with microvessel density suggesting a proangiogenic role in vivo. Thus, MAGP2 may serve as a survival-associated target.}
\usage{data(GSE18520_eset)}
\format{

Authors: Mok SC, Bonome T, Vathipadiekal V, Bell A, Johnson ME, Wong KK, Park DC, Hao K, Yip DK, Donninger H, Ozbun L, Samimi G, Brady J, Randonovich M, Pise-Masison CA, Barrett JC, Wong WH, Welch WR, Berkowitz RS, Birrer MJ.

Lab: Department of Gynecologic Oncology, The University of Texas M.D. Anderson Cancer Center, Houston, TX 77030, USA.

Contact:mbirrer@partners.org

Title: A gene signature predictive for outcome in advanced ovarian cancer identifies a survival factor: microfibril-associated glycoprotein 2.

Url: http://pubget.com/search?q=19962670

PubMedID: 19962670

}
\details{
Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.

frma normalization used from the frma bioconductor package.

assayData: 18708 features, 63 samples.

GEO_platform_accession: GSE18520
platform_summary: GPL570
biomart_ID: HG-U133_Plus_2 affy_hg_u133_plus_2

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

10 observations deleted due to missingness 

records = 53.00  n.max = 53.00  n.start = 53.00   events = 41.00  median = 2.05    0.95LCL = 1.48 0.95UCL = 3.70 

--------------------------------------------------------------------------------


Available sample meta-data: 

--------------------------------------------------------------------------------

alt_sample_name: 
Min. = 312.0   1st Qu. = 395.0  Median = 694.0  Mean = 893.3  3rd Qu. = 1040.0  Max. = 2237.0 

sample_type: 
healthy = 10. tumor = 53.

subtype: 
ser = 53. NA's = 10.

primarysite: 
ov = 63.

summarygrade:
high = 53. NA's = 10.

summarystage: 
late = 53. NA's = 10.

tumorstage:
4 = 53. NA's = 10.

G: 
4 = 53. NA's = 10.

days_to_death: 
Min. = 150 1st Qu. = 450  Median = 630  Mean = 1212 3rd Qu. = 1440  Max. = 4500 NA's = 10 

vital_status: 
deceased = 41. living = 12. NA's = 10.

debulking: 
unknown = 63.

batch: 

2004-03-12 2004-04-08 2004-04-09 2004-07-20 2004-08-12 2004-08-13 2004-09-30 

        20          6          9         11         10          1          6 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE18520
}
\examples{
data(GSE18520_eset)
## maybe str(GSE18520_eset) ; plot(GSE18520_eset) ...
if(require(affy)){
  summary(GSE18520_eset$vital_status)
  time <- GSE18520_eset$days_to_death / 365
  cens <- ifelse(GSE18520_eset$vital_status=="deceased",1,0)
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
