\name{ GSE30161_eset }
\alias{ GSE30161_eset }
\docType{data}
\title{ Multi-gene expression predictors of single drug responses to adjuvant chemotherapy in ovarian carcinoma: predicting platinum resistance. }
\description{ Despite advances in radical surgery and chemotherapy delivery, ovarian cancer is the most lethal gynecologic malignancy. Standard therapy includes treatment with platinum-based combination chemotherapies yet there is no biomarker model to predict their responses to these agents. We here have developed and independently tested our multi-gene molecular predictors for forecasting patients' responses to individual drugs on a cohort of 55 ovarian cancer patients. To independently validate these molecular predictors, we performed microarray profiling on FFPE tumor samples of 55 ovarian cancer patients (UVA-55) treated with platinum-based adjuvant chemotherapy. Genome-wide chemosensitivity biomarkers were initially discovered from the in vitro drug activities and genomic expression data for carboplatin and paclitaxel, respectively. Multivariate predictors were trained with the cell line data and then evaluated with a historical patient cohort. For the UVA-55 cohort, the carboplatin, taxol, and combination predictors significantly stratified responder patients and non-responder patients (p = 0.019, 0.04, 0.014) with sensitivity = 91\%, 96\%, 93 and NPV = 57\%, 67\%, 67\% in pathologic clinical response. The combination predictor also demonstrated a significant survival difference between predicted responders and non-responders with a median survival of 55.4 months vs. 32.1 months. Thus, COXEN single- and combination-drug predictors successfully stratified platinum resistance and taxane response in an independent cohort of ovarian cancer patients based on their FFPE tumor samples. }
\usage{data( GSE30161_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Ferriss JS, Kim Y, Duska L, Birrer M, Levine DA, Moskaluk C,Theodorescu D, Lee JK 
  Laboratory: Ferriss, Lee 2012 
  Contact information:  
  Title: Multi-gene expression predictors of single drug responses to adjuvant chemotherapy in ovarian carcinoma: predicting platinum resistance. 
  URL:  
  PMIDs: 22348014 

  Abstract: A 215 word abstract is available. Use 'abstract' method.
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
   platform_technology:     
      in situ oligonucleotide

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19816 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19816 features, 58 samples
Platform type: hgu133plus2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  58.00   58.00   58.00   36.00    4.19    2.70    6.17 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       58 character character 

sample_type: 
tumor 
   58 

histological_type: 
   Length     Class      Mode 
       58 character character 

summarygrade: 
high  low NA's 
  33   21    4 

summarystage: 
late 
  58 

tumorstage: 
 3  4 
53  5 

substage: 
 a  b  c 
 9 11 38 

grade: 
   1    2    3 NA's 
   2   19   33    4 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  38.00   53.50   62.00   62.57   72.00   85.00 

pltx: 
 y 
58 

tax: 
 n  y 
 4 54 

neo: 
 n 
58 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   12.0   255.2   386.0   742.1   768.2  4208.0 

recurrence_status: 
norecurrence   recurrence         NA's 
           6           48            4 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   49.0   585.2  1010.0  1375.0  2131.0  4208.0 

vital_status: 
deceased   living 
      36       22 

debulking: 
   optimal suboptimal       NA's 
        26         30          2 

batch: 
2009-10-07 2009-10-08 2009-10-09 2009-10-20 
        28         18          8          4 

uncurated_author_metadata: 
   Length     Class      Mode 
       58 character character 

}}
\keyword{datasets}
