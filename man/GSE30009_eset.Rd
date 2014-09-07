\name{ GSE30009_eset }
\alias{ GSE30009_eset }
\docType{data}
\title{ Multidrug resistance-linked gene signature predicts overall survival of patients with primary ovarian serous carcinoma. }
\description{ This study assesses the ability of multidrug resistance (MDR)-associated gene expression patterns to predict survival in patients with newly diagnosed carcinoma of the ovary. The scope of this research differs substantially from that of previous reports, as a very large set of genes was evaluated whose expression has been shown to affect response to chemotherapy.We applied a customized TaqMan low density array, a highly sensitive and specific assay, to study the expression profiles of 380 MDR-linked genes in 80 tumor specimens collected at initial surgery to debulk primary serous carcinoma. The RNA expression profiles of these drug resistance genes were correlated with clinical outcomes.Leave-one-out cross-validation was used to estimate the ability of MDR gene expression to predict survival. Although gene expression alone does not predict overall survival (OS; P = 0.06), four covariates (age, stage, CA125 level, and surgical debulking) do (P = 0.03). When gene expression was added to the covariates, we found an 11-gene signature that provides a major improvement in OS prediction (log-rank statistic P < 0.003). The predictive power of this 11-gene signature was confirmed by dividing high- and low-risk patient groups, as defined by their clinical covariates, into four specific risk groups on the basis of expression levels.This study reveals an 11-gene signature that allows a more precise prognosis for patients with serous cancer of the ovary treated with carboplatin- and paclitaxel-based therapy. These 11 new targets offer opportunities for new therapies to improve clinical outcome in ovarian cancer. }
\usage{data( GSE30009_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Gillet JP, Calcagno AM, Varma S, Davidson B et al. Multidrug resistance-linked gene signature predicts overall survival of patients with primary ovarian serous carcinoma. Clin Cancer Res 2012 Jun 1;18(11):3197-206. 
  Laboratory: Gillet, Gottesman 2012 
  Contact information:  
  Title: Multidrug resistance-linked gene signature predicts overall survival of patients with primary ovarian serous carcinoma. 
  URL:  
  PMIDs: 22492981 

  Abstract: A 244 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      TaqMan qRT-PCR Homo sapiens Low-Density Array 380
   platform_shorttitle:     
      TaqMan qRT-PCR
   platform_summary:     
      NA
   platform_manufacturer:     
      TaqMan
   platform_distribution:     
      custom
   platform_accession:     
      GPL13728
   platform_technology:     
      qRT-PCR

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ABCA1 ABCA10 ... XRCC6 (359 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 359 features, 103 samples
Platform type: NA
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 103.00  103.00  103.00   57.00    3.42    2.92    5.34 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      103 character character 

sample_type: 
borderline      tumor 
         1        102 

histological_type: 
clearcell       ser 
        1       102 

summarygrade: 
high  low NA's 
  92    9    2 

summarystage: 
late 
 103 

tumorstage: 
 3  4 
82 21 

substage: 
   b    c NA's 
   2   60   41 

grade: 
   1    2    3 NA's 
   4    5   92    2 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  30.00   56.00   61.00   62.45   71.50   87.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     24     598    1053    1156    1568    4748 

vital_status: 
deceased   living 
      57       46 

debulking: 
   optimal suboptimal 
        81         22 

uncurated_author_metadata: 
   Length     Class      Mode 
      103 character character 

}}
\keyword{datasets}
