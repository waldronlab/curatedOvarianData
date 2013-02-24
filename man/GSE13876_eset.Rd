\name{ GSE13876_eset }
\alias{ GSE13876_eset }
\docType{data}
\title{ Survival-related profile, pathways, and transcription factors in ovarian cancer. }
\description{ Ovarian cancer has a poor prognosis due to advanced stage at presentation and either intrinsic or acquired resistance to classic cytotoxic drugs such as platinum and taxoids. Recent large clinical trials with different combinations and sequences of classic cytotoxic drugs indicate that further significant improvement in prognosis by this type of drugs is not to be expected. Currently a large number of drugs, targeting dysregulated molecular pathways in cancer cells have been developed and are introduced in the clinic. A major challenge is to identify those patients who will benefit from drugs targeting these specific dysregulated pathways.The aims of our study were (1) to develop a gene expression profile associated with overall survival in advanced stage serous ovarian cancer, (2) to assess the association of pathways and transcription factors with overall survival, and (3) to validate our identified profile and pathways/transcription factors in an independent set of ovarian cancers.According to a randomized design, profiling of 157 advanced stage serous ovarian cancers was performed in duplicate using approximately 35,000 70-mer oligonucleotide microarrays. A continuous predictor of overall survival was built taking into account well-known issues in microarray analysis, such as multiple testing and overfitting. A functional class scoring analysis was utilized to assess pathways/transcription factors for their association with overall survival. The prognostic value of genes that constitute our overall survival profile was validated on a fully independent, publicly available dataset of 118 well-defined primary serous ovarian cancers. Furthermore, functional class scoring analysis was also performed on this independent dataset to assess the similarities with results from our own dataset. An 86-gene overall survival profile discriminated between patients with unfavorable and favorable prognosis (median survival, 19 versus 41 mo, respectively; permutation p-value of log-rank statistic = 0.015) and maintained its independent prognostic value in multivariate analysis. Genes that composed the overall survival profile were also able to discriminate between the two risk groups in the independent dataset. In our dataset 17/167 pathways and 13/111 transcription factors were associated with overall survival, of which 16 and 12, respectively, were confirmed in the independent dataset.Our study provides new clues to genes, pathways, and transcription factors that contribute to the clinical outcome of serous ovarian cancer and might be exploited in designing new treatment strategies. }
\usage{data( GSE13876_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Crijns AP, Fehrmann RS, de Jong S, Gerbens F, Meersma GJ, Klip HG, Hollema H, Hofstra RM, te Meerman GJ, de Vries EG, van der Zee AG.Survival-related profile, pathways, and transcription factors in ovarian cancer. PLoS Med. 2009 Feb 3; 6(2):e24. 
  Laboratory: Crijns, van der Zee 2009 
  Contact information:  
  Title: Survival-related profile, pathways, and transcription factors in ovarian cancer. 
  URL:  
  PMIDs: 19192944 

  Abstract: A 371 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Operon human v3 ~35K 70-mer two-color oligonucleotide microarrays
   platform_shorttitle:     
      Operon v3 two-color
   platform_summary:     
      OperonHumanV3
   platform_manufacturer:     
      other
   platform_distribution:     
      non-commercial
   platform_accession:     
      GPL7759

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (20505 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20505 features, 157 samples
Platform type: OperonHumanV3
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 157.00  157.00  157.00  113.00    2.05    1.56    2.71 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
 151 NA's 
   1  156 

unique_patient_ID: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      1      40      79      79     118     157 

sample_type: 
tumor 
  157 

histological_type: 
ser 
157 

primarysite: 
 ov 
157 

summarygrade: 
high  low NA's 
  85   59   13 

summarystage: 
late 
 157 

grade: 
   1    2    3    4 NA's 
  14   45   82    3   13 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   50.00   60.00   57.95   67.00   84.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30     360     630    1100    1470    7020 

vital_status: 
deceased   living 
     113       44 

uncurated_author_metadata: 
   Length     Class      Mode 
      157 character character 

}}
\keyword{datasets}
