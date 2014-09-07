\name{ GSE32063_eset }
\alias{ GSE32063_eset }
\docType{data}
\title{ High-risk ovarian cancer based on 126-gene expression signature is uniquely characterized by downregulation of antigen presentation pathway. }
\description{ High-grade serous ovarian cancers are heterogeneous not only in terms of clinical outcome but also at the molecular level. Our aim was to establish a novel risk classification system based on a gene expression signature for predicting overall survival, leading to suggesting novel therapeutic strategies for high-risk patients.In this large-scale cross-platform study of six microarray data sets consisting of 1,054 ovarian cancer patients, we developed a gene expression signature for predicting overall survival by applying elastic net and 10-fold cross-validation to a Japanese data set A (n = 260) and evaluated the signature in five other data sets. Subsequently, we investigated differences in the biological characteristics between high- and low-risk ovarian cancer groups.An elastic net analysis identified a 126-gene expression signature for predicting overall survival in patients with ovarian cancer using the Japanese data set A (multivariate analysis, P = 4 ?? 10(-20)). We validated its predictive ability with five other data sets using multivariate analysis (Tothill's data set, P = 1 ?? 10(-5); Bonome's data set, P = 0.0033; Dressman's data set, P = 0.0016; TCGA data set, P = 0.0027; Japanese data set B, P = 0.021). Through gene ontology and pathway analyses, we identified a significant reduction in expression of immune-response-related genes, especially on the antigen presentation pathway, in high-risk ovarian cancer patients.This risk classification based on the 126-gene expression signature is an accurate predictor of clinical outcome in patients with advanced stage high-grade serous ovarian cancer and has the potential to develop new therapeutic strategies for high-grade serous ovarian cancer patients. }
\usage{data( GSE32063_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Yoshihara K, Tsunoda T, Shigemizu D, Fujiwara H et al. High-risk ovarian cancer based on 126-gene expression signature is uniquely characterized by downregulation of antigen presentation pathway. Clin Cancer Res 2012 Mar 1;18(5):1374-85. 
  Laboratory: Yoshihara, Tanaka 2012 
  Contact information:  
  Title: High-risk ovarian cancer based on 126-gene expression signature is uniquely characterized by downregulation of antigen presentation pathway. 
  URL:  
  PMIDs: 22241791 

  Abstract: A 255 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name vers
ion)
   platform_shorttitle:     
      Agilent G4112F
   platform_summary:     
      hgug4112a
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      GPL6480
   platform_technology:     
      in situ oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (20106 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20106 features, 40 samples
Platform type: hgug4112a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  40.00   40.00   40.00   22.00    4.44    3.29      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       40 character character 

sample_type: 
tumor 
   40 

histological_type: 
ser 
 40 

summarygrade: 
high  low 
  17   23 

summarystage: 
late 
  40 

tumorstage: 
 3  4 
31  9 

substage: 
   b    c NA's 
   3   28    9 

grade: 
 2  3 
23 17 

pltx: 
 y 
40 

tax: 
 y 
40 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    210     705    1155    1346    1792    3330 

vital_status: 
deceased   living 
      22       18 

debulking: 
   optimal suboptimal 
        19         21 

uncurated_author_metadata: 
   Length     Class      Mode 
       40 character character 

}}
\keyword{datasets}
