\name{ PMID19318476_eset }
\alias{ PMID19318476_eset }
\docType{data}
\title{ Microarray analysis of early stage serous ovarian cancers shows profiles predictive of favorable outcome. }
\description{ Although few women with advanced serous ovarian cancer are cured, detection of the disease at an early stage is associated with a much higher likelihood of survival. We previously used gene expression array analysis to distinguish subsets of advanced cancers based on disease outcome. In the present study, we report on gene expression of early-stage cancers and validate our prognostic model for advanced-stage cancers.Frozen specimens from 39 stage I/II, 42 stage III/IV, and 20 low malignant potential cancers were obtained from four different sites. A linear discriminant model was used to predict survival based upon array data.We validated the late-stage survival model and show that three of the most differentially expressed genes continue to be predictive of outcome. Most early-stage cancers (38 of 39 invasive, 15 of 20 low malignant potential) were classified as long-term survivors (median probabilities 0.97 and 0.86). MAL, the most differentially expressed gene, was further validated at the protein level and found to be an independent predictor of poor survival in an unselected group of advanced serous cancers (P = 0.0004).These data suggest that serous ovarian cancers detected at an early stage generally have a favorable underlying biology similar to advanced-stage cases that are long-term survivors. Conversely, most late-stage ovarian cancers seem to have a more virulent biology. This insight suggests that if screening approaches are to succeed it will be necessary to develop approaches that are able to detect these virulent cancers at an early stage. }
\usage{data( PMID19318476_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Berchuck A, Iversen ES, Luo J, Clarke JP, Horne H, Levine DA, Boyd J, Alonso MA, Secord AA, Bernardini MQ, Barnett JC, Boren T, Murphy SK, Dressman HK, Marks JR, Lancaster JM.Microarray analysis of early stage serous ovarian cancers shows profiles predictive of favorable outcome. Clin Cancer Res. 2009 Apr 1; 15(7):2448-55. 
  Laboratory: Berchuck, Lancaster 2009 
  Contact information:  
  Title: Microarray analysis of early stage serous ovarian cancers shows profiles predictive of favorable outcome. 
  URL:  
  PMIDs: 19318476 

  Abstract: A 241 word abstract is available. Use 'abstract' method.
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
   warnings:     
      These samples are a subset of PMID17290060.

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12981 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12981 features, 42 samples
Platform type: hgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  42.00   42.00   42.00   22.00    2.79    2.30      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       42 character character 

sample_type: 
tumor 
   42 

histological_type: 
ser 
 42 

summarygrade: 
high  low NA's 
  24   17    1 

summarystage: 
early  late  NA's 
    2    39     1 

tumorstage: 
   1    2    3    4 NA's 
   1    1   29   10    1 

substage: 
   a    b    c NA's 
   1    1   29   11 

grade: 
   1    2    3 NA's 
   2   15   24    1 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  33.00   55.00   62.00   61.46   70.00   81.00       1 

recurrence_status: 
norecurrence   recurrence 
           6           36 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   30.0   367.5   825.0  1105.0  1050.0  3420.0 

vital_status: 
deceased   living 
      22       20 

debulking: 
   optimal suboptimal       NA's 
        20         21          1 

batch: 
2004-03-09 2004-03-16 2004-04-20 2004-05-18 2004-05-21 2004-05-27 2004-06-22 
        14          3          4          8          6          5          1 
2004-06-23 
         1 

uncurated_author_metadata: 
   Length     Class      Mode 
       42 character character 

}}
\keyword{datasets}
