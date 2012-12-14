\name{GSE17260_eset}
\alias{GSE17260_eset}
\docType{data}
\title{Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets.}
\description{Advanced-stage ovarian cancer patients are generally treated with platinum/taxane-based chemotherapy after primary debulking surgery. However, there is a wide range of outcomes for individual patients. Therefore, the clinicopathological factors alone are insufficient for predicting prognosis. Our aim is to identify a progression-free survival (PFS)-related molecular profile for predicting survival of patients with advanced-stage serous ovarian cancer.}
\usage{data(GSE17260_eset)}
\format{

Authors: Yoshihara K, Tajima A, Yahata T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo Y, Kotera K, Masuzaki H, Tashiro H, Katabuchi H, Inoue I, Tanaka K.

Lab : Department of Obstetrics and Gynecology, Niigata University Graduate School of Medical and Dental Sciences, Niigata, Japan.

Contact: tanaken@med.niigata-u.ac.jp

Title: Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets.

Url: http://pubget.com/search?q=20300634

PubMedID: 20300634
}

\details{
Platform used: Agilent-012391 Whole Human Genome Oligo Microarray G4112A.

Data normalization was performed in GeneSpring GX 10 (Agilent Technologies).  Processed data here are those provided by the original study authors.

assayData: 18844 features, 110 samples.

GEO_platform_accession: GSE17260.
platform_summary = GPL6848.
biomart_ID = Agilent G4112A agilent_wholegenome.

Overall survival time-to-event summary (in years): 
Call: survfit(formula = Surv(time, cens) ~ -1)

records = 110.00  n.max = 110.00  n.start = 110.00   events = 46.00   median = 4.44 0.95LCL = 4.03 0.95UCL = NA

---------------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------------

alt_sample_name:
Length = 110.
Class = character.
Mode = character.

sample_type: 
tumor = 110.

subtype: 
ser = 110.

primarysite: 
ov = 110.

summarystage: 
late = 110.

tumorstage:
3 = 93. 4 = 17.

substage: 
a = 6. b = 18. c = 69. NA's = 17.

pltx: 
y = 110.

tax:
y = 110.

days_to_tumor_recurrence: 
Min. = 30.0  1st Qu. = 285.0  Median = 510.0   Mean = 673.9  3rd Qu. = 870.0  Max. = 2250.0 

recurrence_status: 
norecurrence = 34.   recurrence = 76.

days_to_death: 
Min. = 30 1st Qu. =  660 Median = 915 Mean = 1086 3rd Qu. = 1530 Max. = 2430

vital_status: 
deceased = 46. living = 64.

debulking: 
optimal = 57. suboptimal = 53.
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17260
}
\examples{
data(GSE17260_eset)
## maybe str(GSE17260_eset) ; plot(GSE17260_eset) ...
if(require(affy)){
  summary(GSE17260_eset$recurrence_status)
  time <- GSE17260_eset$days_to_death / 365
  cens <- ifelse(GSE17260_eset$vital_status=="deceased",1,0)
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
