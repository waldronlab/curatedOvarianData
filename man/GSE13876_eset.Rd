\name{GSE13876_eset}
\alias{GSE13876_eset}
\docType{data}
\title{Survival-related profile, pathways, and transcription factors in ovarian cancer.
}
\description{Ovarian cancer has a poor prognosis due to advanced stage at presentation and either intrinsic or acquired resistance to classic cytotoxic drugs such as platinum and taxoids. Recent large clinical trials with different combinations and sequences of classic cytotoxic drugs indicate that further significant improvement in prognosis by this type of drugs is not to be expected. Currently a large number of drugs, targeting dysregulated molecular pathways in cancer cells have been developed and are introduced in the clinic. A major challenge is to identify those patients who will benefit from drugs targeting these specific dysregulated pathways.The aims of our study were (1) to develop a gene expression profile associated with overall survival in advanced stage serous ovarian cancer, (2) to assess the association of pathways and transcription factors with overall survival, and (3) to validate our identified profile and pathways/transcription factors in an independent set of ovarian cancers.}
\usage{data(GSE13876_eset)}
\format{

Authors: Crijns AP, Fehrmann RS, de Jong S, Gerbens F, Meersma GJ, Klip HG, Hollema H, Hofstra RM, te Meerman GJ, de Vries EG, van der Zee AG.

Lab: Department of Gynecologic Oncology, University Medical Centre Groningen, University of Groningen, Groningen, The Netherlands.

Contact: A.G.J.van.der.Zee@og.umcg.nl

Title: Survival-Related Profile, Pathways, and Transcription Factors in Ovarian Cancer.

Url: http://pubget.com/search?q=19192944

PubMedID: 19192944
}
\details{
Platform used: Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays.

Quantile normalization was used.  Processed data here are those provided by the original study authors.

assayData: 20112 features, 157 samples.

GEO_platform_accession: GSE13876.
platform_summary : GPL7759.
biomart_ID: Operon two-color oligo       

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)


records: 157.00   n.max = 157.00 n.start = 157.00  events = 113.00 median = 2.05 0.95LCL = 1.56 0.95UCL = 2.71
          
--------------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------------
alt_sample_name: 
Length = 157.
Class = character.      
Mode = character.

unique_patient_ID: 
Min. = 1. 1st Qu. = 40. Median = 79. Mean = 79. 3rd Qu. = 118. Max. = 157.

sample_type: 
tumor = 157.

subtype:
ser = 157.
157 

primarysite: 
ov = 157.

summarystage: 
late = 157.
 
tumorstage
4 = 157.

age_at_initial_pathologic_diagnosis: 
Min. = 21.00 1st Qu. = 50.00 Median = 60.00 Mean = 57.95 3rd Qu. = 67.00 Max = 84.00 

recurrence_status: 
norecurrence = 44.  recurrence = 113.

days_to_death: 
Min. = 30. 1st Qu. = 360. Median = 630. Mean = 1100. 3rd Qu. = 1470. Max. = 7020.

vital_status:
deceased = 113. living = 44.

debulking: 
unknown = 157.
    157 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE13876
}
\examples{
data(GSE13876_eset)
## maybe str(GSE13876_eset) ; plot(GSE13876_eset) ...
if(require(affy)){
  summary(GSE13876_eset$recurrence_status)
  time <- GSE13876_eset$days_to_death / 365
  cens <- ifelse(GSE13876_eset$vital_status=="deceased",1,0)
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
