\name{TCGA.mirna.8x15kv2_eset}
\alias{TCGA.mirna.8x15kv2_eset}
\docType{data}
\title{Ovarian serous cystadenocarcinoma: Sample Counts and Findings}
\description{The Cancer Genome Atlas}
\usage{data(TCGA.mirna.8x15kv2_eset)}
\format{

Authors: The Cancer Genome Atlas.

}
\details{
\preformatted{

Level 3 (fully-processed) microRNA data downloaded from TCGA.

assayData: 799 features, 554 samples

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

   10 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 544.00  544.00  544.00  286.00    3.71    3.42    4.03 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      554 character character 

sample_type: 
tumor 
  554 

subtype: 
ser 
554 

primarysite: 
other    ov 
    4   550 

summarygrade: 
high  low NA's 
 474   68   12 

summarystage: 
early  late  NA's 
   39   511     4 

tumorstage: 
   1    2    3    4 NA's 
  16   23  426   85    4 

substage: 
   b    c NA's 
  31  434   89 

G: 
   1    2    3    4 NA's 
   6   62  473    1   12 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   51.00   59.00   59.81   69.00   89.00 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   235.2   436.0   618.7   797.0  5480.0      44 

recurrence_status: 
norecurrence   recurrence 
         268          286 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   346.0   867.5   997.7  1446.0  5480.0      10 

vital_status: 
deceased   living     NA's 
     286      261        7 

site_of_tumor_first_recurrence: 
                locoregional locoregional_plus_metastatic 
                         145                            3 
                  metastasis                         NA's 
                         138                          268 

primary_therapy_outcome_success: 
  completeresponse    partialresponse progressivedisease      stabledisease 
               308                 63                 41                 30 
              NA's 
               112 

debulking: 
   optimal suboptimal       NA's 
       112        384         58 
}


}
\examples{
data(TCGA.mirna.8x15kv2_eset)
}
\keyword{datasets}
