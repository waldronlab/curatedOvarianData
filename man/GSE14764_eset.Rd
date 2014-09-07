\name{ GSE14764_eset }
\alias{ GSE14764_eset }
\docType{data}
\title{ A prognostic gene expression index in ovarian cancer - validation across different independent data sets. }
\description{ Ovarian carcinoma has the highest mortality rate among gynaecological malignancies. In this project, we investigated the hypothesis that molecular markers are able to predict outcome of ovarian cancer independently of classical clinical predictors, and that these molecular markers can be validated using independent data sets. We applied a semi-supervised method for prediction of patient survival. Microarrays from a cohort of 80 ovarian carcinomas (TOC cohort) were used for the development of a predictive model, which was then evaluated in an entirely independent cohort of 118 carcinomas (Duke cohort). A 300-gene ovarian prognostic index (OPI) was generated and validated in a leave-one-out approach in the TOC cohort (Kaplan-Meier analysis, p = 0.0087). In a second validation step, the prognostic power of the OPI was confirmed in an independent data set (Duke cohort, p = 0.0063). In multivariate analysis, the OPI was independent of the post-operative residual tumour, the main clinico-pathological prognostic parameter with an adjusted hazard ratio of 6.4 (TOC cohort, CI 1.8-23.5, p = 0.0049) and 1.9 (Duke cohort, CI 1.2-3.0, p = 0.0068). We constructed a combined score of molecular data (OPI) and clinical parameters (residual tumour), which was able to define patient groups with highly significant differences in survival. The integrated analysis of gene expression data as well as residual tumour can be used for optimized assessment of the prognosis of platinum-taxol-treated ovarian cancer. As traditional treatment options are limited, this analysis may be able to optimize clinical management and to identify those patients who would be candidates for new therapeutic strategies. }
\usage{data( GSE14764_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Denkert C, Budczies J, Darb-Esfahani S, Gy??rffy B et al. A prognostic gene expression index in ovarian cancer - validation across different independent data sets. J Pathol 2009 Jun;218(2):273-80. 
  Laboratory: Denkert, Lage 2009 
  Contact information:  
  Title: A prognostic gene expression index in ovarian cancer - validation across different independent data sets. 
  URL:  
  PMIDs: 19294737 

  Abstract: A 254 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133A] Affymetrix Human Genome U133A Array
   platform_shorttitle:     
      Affymetrix HG-U133A
   platform_summary:     
      hgu133a
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL96
   platform_technology:     
      in situ oligonucleotide

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (13104 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 13104 features, 80 samples
Platform type: hgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  80.00   80.00   80.00   21.00    4.52    4.19      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1.00   20.75   40.50   40.50   60.25   80.00 

sample_type: 
tumor 
   80 

histological_type: 
   Length     Class      Mode 
       80 character character 

primarysite: 
ov 
80 

summarygrade: 
high  low 
  54   26 

summarystage: 
early  late 
    9    71 

tumorstage: 
 1  2  3  4 
 8  1 69  2 

substage: 
   a    b    c NA's 
   4    6   32   38 

grade: 
 1  2  3 
 3 23 54 

recurrence_status: 
norecurrence   recurrence         NA's 
          50           26            4 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    210     660    1050    1011    1328    2190 

vital_status: 
deceased   living 
      21       59 

batch: 
   Length     Class      Mode 
       80 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
       80 character character 

}}
\keyword{datasets}
