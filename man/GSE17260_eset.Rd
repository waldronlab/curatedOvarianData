\name{ GSE17260_eset }
\alias{ GSE17260_eset }
\docType{data}
\title{ Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets. }
\description{ Advanced-stage ovarian cancer patients are generally treated with platinum/taxane-based chemotherapy after primary debulking surgery. However, there is a wide range of outcomes for individual patients. Therefore, the clinicopathological factors alone are insufficient for predicting prognosis. Our aim is to identify a progression-free survival (PFS)-related molecular profile for predicting survival of patients with advanced-stage serous ovarian cancer.Advanced-stage serous ovarian cancer tissues from 110 Japanese patients who underwent primary surgery and platinum/taxane-based chemotherapy were profiled using oligonucleotide microarrays. We selected 88 PFS-related genes by a univariate Cox model (p<0.01) and generated the prognostic index based on 88 PFS-related genes after adjustment of regression coefficients of the respective genes by ridge regression Cox model using 10-fold cross-validation. The prognostic index was independently associated with PFS time compared to other clinical factors in multivariate analysis [hazard ratio (HR), 3.72; 95\% confidence interval (CI), 2.66-5.43; p<0.0001]. In an external dataset, multivariate analysis revealed that this prognostic index was significantly correlated with PFS time (HR, 1.54; 95\% CI, 1.20-1.98; p = 0.0008). Furthermore, the correlation between the prognostic index and overall survival time was confirmed in the two independent external datasets (log rank test, p = 0.0010 and 0.0008).The prognostic ability of our index based on the 88-gene expression profile in ridge regression Cox hazard model was shown to be independent of other clinical factors in predicting cancer prognosis across two distinct datasets. Further study will be necessary to improve predictive accuracy of the prognostic index toward clinical application for evaluation of the risk of recurrence in patients with advanced-stage serous ovarian cancer. }
\usage{data( GSE17260_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Yoshihara K, Tajima A, Yahata T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo Y, Kotera K, Masuzaki H, Tashiro H, Katabuchi H, Inoue I, Tanaka K.Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets. PLoS One. 2010 Mar 12; 5(3):e9615. 
  Laboratory: Yoshihara, Tanaka 2010 
  Contact information:  
  Title: Gene expression profile for predicting survival in advanced-stage serous ovarian cancer across two independent datasets. 
  URL:  
  PMIDs: 20300634 

  Abstract: A 257 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent-012391 Whole Human Genome Oligo Microarray G4112A
   platform_shorttitle:     
      Agilent G4112A
   platform_summary:     
      hgug4112a
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      GPL6848
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
assayData: 20106 features, 110 samples
Platform type: hgug4112a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 110.00  110.00  110.00   46.00    4.44    4.03      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      110 character character 

sample_type: 
tumor 
  110 

histological_type: 
ser 
110 

primarysite: 
 ov 
110 

summarygrade: 
high  low 
  43   67 

summarystage: 
late 
 110 

tumorstage: 
 3  4 
93 17 

substage: 
   a    b    c NA's 
   6   18   69   17 

grade: 
 1  2  3 
26 41 43 

pltx: 
  y 
110 

tax: 
  y 
110 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   30.0   285.0   510.0   673.9   870.0  2250.0 

recurrence_status: 
norecurrence   recurrence 
          34           76 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30     660     915    1086    1530    2430 

vital_status: 
deceased   living 
      46       64 

debulking: 
   optimal suboptimal 
        57         53 

uncurated_author_metadata: 
   Length     Class      Mode 
      110 character character 

}}
\keyword{datasets}
