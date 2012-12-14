\name{TCGA_eset}
\alias{TCGA_eset}
\docType{data}
\title{Ovarian serous cystadenocarcinoma: Sample Counts and Findings}
\description{The Cancer Genome Atlas}
\usage{data(TCGA_eset)}
\format{

Authors: The Cancer Genome Atlas.

Title: Ovarian serous cystadenocarcinoma: Sample Counts and Findings.

Url: http://tcga-data.nci.nih.gov/tcga/tcgaCancerDetails.jsp?diseaseType=OV&diseaseName=Ovarian\%20serous\%20cystadenocarcinoma

}
\details{

  Note that the TCGA dataset contains some samples replicated across
  batches, but there are no clinical metadata associated with these
  samples.  The following provides these sample pairs and the Pearson
  correlation of their RMA-normalized expression profiles:

  \preformatted{
  name1         name2	      pearson correlation 
  TCGA.24.2290  TCGA.13.2071  0.9969045 
  TCGA.25.2392  TCGA.29.2432  0.9976206 
  TCGA.25.2404  TCGA.36.2537  0.9978256 
  TCGA.24.2262  TCGA.42.2591  0.9967085 
  TCGA.24.2293  TCGA.23.2641  0.997569 
  TCGA.25.2393  TCGA.29.2434  0.996287 
  TCGA.59.2350  TCGA.36.2548  0.9968095 
  TCGA.24.2295  TCGA.23.2643  0.9952585 
  TCGA.09.2048  TCGA.04.1644  0.9967052 
  TCGA.24.2271  TCGA.13.2057  0.9969408 
  TCGA.25.2397  TCGA.61.2612  0.9957279 
  TCGA.09.2049  TCGA.42.2582  0.9932573 
  TCGA.24.2280  TCGA.13.2059  0.9937523 
  TCGA.24.2298  TCGA.23.2647  0.9970841 
  TCGA.25.2398  TCGA.36.2530  0.996212 
  TCGA.09.2051  TCGA.42.2588  0.9971348 
  TCGA.24.2288  TCGA.13.2065  0.9954178 
  TCGA.29.2428  TCGA.36.2544  0.9964141 
  TCGA.24.2289  TCGA.13.2066  0.9968514 
  TCGA.36.2529  TCGA.29.1704  0.9976284 
  TCGA.36.2532  TCGA.24.1852  0.997482 
}

rma normalization used from the rma bioconductor package.

assayData: 12858 features, 570 samples.

GEO_platform_accession: TCGA

platform_summary: GPL4685

biomart_ID: HT-HG_U133A affy_hg_u133a

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

records = 576.  n.max = 576. n.start = 576.  events = 0. median = NA. 0.95LCL = NA. 0.95UCL = NA.

-------------------------------------------------------------------------------


Available sample meta-data: 

-------------------------------------------------------------------------------

alt_sample_name: 
Length = 576.
Class = character.
Mode = character.

unique_patient_ID: 
Length = 576.
Class = character.
Mode = character.

sample_type: 
tumor = 576.

subtype: 
ser = 576.

primarysite: 
ov = 576.

G: 
1 = 6.
2 = 68.
3 = 485.
4 = 1.

days_to_death: 
Min. = 8  Median = 1024   Mean = 1083.22  Max. = 4623

batch: 
Length = 570.
Class = character.
Mode = character.

}
\examples{
data(TCGA_eset)
## maybe str(TCGA_eset) ; plot(TCGA_eset) ...

if(require(affy)){
  summary(TCGA_eset$recurrence_status)
  time <- TCGA_eset$days_to_death / 365
  cens <- ifelse(TCGA_eset$vital_status=="deceased",1,0)
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
