\name{ GSE9891_eset }
\alias{ GSE9891_eset }
\docType{data}
\title{ Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome. }
\description{ The study aim to identify novel molecular subtypes of ovarian cancer by gene expression profiling with linkage to clinical and pathologic features.Microarray gene expression profiling was done on 285 serous and endometrioid tumors of the ovary, peritoneum, and fallopian tube. K-means clustering was applied to identify robust molecular subtypes. Statistical analysis identified differentially expressed genes, pathways, and gene ontologies. Laser capture microdissection, pathology review, and immunohistochemistry validated the array-based findings. Patient survival within k-means groups was evaluated using Cox proportional hazards models. Class prediction validated k-means groups in an independent dataset. A semisupervised survival analysis of the array data was used to compare against unsupervised clustering results.Optimal clustering of array data identified six molecular subtypes. Two subtypes represented predominantly serous low malignant potential and low-grade endometrioid subtypes, respectively. The remaining four subtypes represented higher grade and advanced stage cancers of serous and endometrioid morphology. A novel subtype of high-grade serous cancers reflected a mesenchymal cell type, characterized by overexpression of N-cadherin and P-cadherin and low expression of differentiation markers, including CA125 and MUC1. A poor prognosis subtype was defined by a reactive stroma gene expression signature, correlating with extensive desmoplasia in such samples. A similar poor prognosis signature could be found using a semisupervised analysis. Each subtype displayed distinct levels and patterns of immune cell infiltration. Class prediction identified similar subtypes in an independent ovarian dataset with similar prognostic trends.Gene expression profiling identified molecular subtypes of ovarian cancer of biological and clinical importance. }
\usage{data( GSE9891_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Tothill RW, Tinker AV, George J, Brown R, Fox SB, Lade S, Johnson DS, Trivett MK, Etemadmoghadam D, Locandro B, Traficante N, Fereday S, Hung JA, Chiew YE, Haviv I, Gertig D, DeFazio A, Bowtell DD.Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome. Clin Cancer Res. 2008 Aug 15; 14(16):5198-208. 
  Laboratory: Tothill, Bowtell 2008 
  Contact information:  
  Title: Novel molecular subtypes of serous and endometrioid ovarian cancer linked to clinical outcome. 
  URL:  
  PMIDs: 18698038 

  Abstract: A 243 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      Affymetrix HG-U133Plus2
   platform_summary:     
      hgu133plus2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL570

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (19093 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19093 features, 285 samples
Platform type: hgu133plus2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   7 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 278.00  278.00  278.00  113.00    3.95    3.53    5.01 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      285 character character 

sample_type: 
tumor 
  285 

histological_type: 
 endo other   ser 
   20     1   264 

primarysite: 
   ft other    ov 
    8    34   243 

arrayedsite: 
   ft other    ov 
    2    83   200 

summarygrade: 
high  low NA's 
 163  116    6 

summarystage: 
early  late  NA's 
   42   240     3 

tumorstage: 
   1    2    3    4 NA's 
  24   18  218   22    3 

substage: 
   a    b    c NA's 
  26   19  212   28 

grade: 
   1    2    3 NA's 
  19   97  163    6 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  22.00   53.00   59.00   59.62   68.00   80.00       3 

pltx: 
   n    y NA's 
  39  243    3 

tax: 
   n    y NA's 
  87  195    3 

neo: 
   n    y NA's 
 264   18    3 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    0.0   300.0   450.0   618.9   810.0  4980.0      10 

recurrence_status: 
norecurrence   recurrence         NA's 
          94          188            3 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    0.0   547.5   855.0   955.1  1252.0  6420.0       7 

vital_status: 
deceased   living     NA's 
     113      169        3 

debulking: 
   optimal suboptimal       NA's 
       160         88         37 

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

uncurated_author_metadata: 
   Length     Class      Mode 
      285 character character 

}}
\keyword{datasets}
