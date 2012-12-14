\name{GSE9891_eset}
\alias{GSE9891_eset}
\docType{data}
\title{Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome.}
\description{The study aim to identify novel molecular subtypes of ovarian cancer by gene expression profiling with linkage to clinical and pathologic features.
}
\usage{data(GSE9891_eset)}
\format{

Authors: Tothill RW, Tinker AV, George J, Brown R, Fox SB, Lade S, Johnson DS, Trivett MK, Etemadmoghadam D, Locandro B, Traficante N, Fereday S, Hung JA, Chiew YE, Haviv I, Gertig D, DeFazio A, Bowtell DD.

Lab: Peter MacCallum Cancer Center, University of Melbourne, Melbourne, Australia.

Contact: d.bowtell@petermac.org

Title: Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome.

Url: http://pubget.com/search?q=18698038

PubMedID: 18698038
}
\details{
Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.

frma normalization used from the frma bioconductor package.

assayData: 18708 features, 285 samples.

GEO_platform_accession: GSE9891

platform_summary: GPL570

biomart_ID: HG-U133_Plus_2 affy_hg_u133_plus_2

Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

10 observations deleted due to missingness 

records = 275.00  n.max = 275.00 n.start = 275.00  events = 113.00 median = 3.95 0.95LCL = 3.53 0.95UCL = 5.01

---------------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------------

alt_sample_name: 
Length = 285.
Class = character.
Mode = character.

sample_type:
tumor = 285.

subtype:
endo = 20. other = 1. ser = 264.
   
primarysite: 
ft = 8. other = 34. ov = 243.

arrayedsite: 
ft = 2. other = 83. ov = 200.
    
summarygrade: 
high = 163. low = 116. NA's = 6.

summarystage: 
early = 42. late = 240. NA's = 3.

tumorstage:
1 = 24. 2 = 18. 3 = 218. 4 = 22. NA's = 3.

substage: 
a = 26.    b = 19. c = 212. NA's = 28.

G:
1 = 19. 2 = 97. 3 = 163. NA's = 6.

age_at_initial_pathologic_diagnosis: 
Min. = 22.00 1st Qu. = 53.00 Median = 59.00  Mean = 59.62 3rd Qu. = 68.00 Max. = 80.00 NA's = 3.00

pltx: 
n = 39. y = 243. NA's = 3.

tax: 
n = 87. y = 195. NA's = 3.

neo:
n = 264.  y = 18. NA's = 3.
 
days_to_tumor_recurrence: 
   Min. = 0.0 1st Qu. = 300.0  Median = 450.0   Mean = 621.6 3rd Qu. = 810.0 Max. = 4980.0 NA's = 7.0

recurrence_status:
norecurrence = 94. recurrence = 188. NA's = 3.

days_to_death: 
   Min. = 0.0 1st Qu. = 555.0  Median = 870.0 Mean = 956.1 3rd Qu. = 1245.0    Max. = 6420.0 NA's = 10.0

vital_status: 
deceased = 113. living = 169. NA's = 3.

debulking: 
optimal = 160. suboptimal = 88.    unknown = 37.

batch: 

2004-12-03 2004-12-23 2005-01-12 2005-01-17 2005-01-24 2005-01-31 2005-02-21 

         3          4          7          7          8         10         10 
         
2005-03-17 2005-05-05 2005-05-09 2005-05-25 2005-05-27 2005-05-30 2005-06-02 

         2          1          1          2          3          3          6 
         
2005-06-06 2005-06-08 2005-06-16 2005-06-17 2005-06-24 2005-07-06 2005-07-15 

         4          5          3          5          6          2          9 
         
2005-07-20 2005-07-29 2005-08-03 2005-08-05 2005-08-18 2005-08-24 2005-08-26 

         7          5          6          3          4          8          4 
         
2005-09-09 2005-09-14 2005-09-16 2005-09-21 2005-10-05 2005-10-26 2005-10-28 

         4          6          6          4          5          2          4 
         
2005-11-04 2005-11-09 2005-11-11 2005-11-23 2005-12-15 2005-12-21 2006-01-20 

         6          3          7          4          7          8          3 
         
2006-01-31 2006-02-08 2006-02-28 2006-04-05 2006-04-06 2006-04-12 2006-04-13 

         7          3          3          7          3          7          4 
         
2006-04-28 2006-05-03 2006-06-06 2006-06-07 2006-06-22 2006-07-07 2006-07-19 

         6          9          6          3          9          4          7 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE9891
}
\examples{
data(GSE9891_eset)
if(require(affy)){
  summary(GSE9891_eset$recurrence_status)
  time <- GSE9891_eset$days_to_death / 365
  cens <- ifelse(GSE9891_eset$vital_status=="deceased",1,0)
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
