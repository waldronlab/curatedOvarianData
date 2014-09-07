\name{ TCGA.RNASeqV2_eset }
\alias{ TCGA.RNASeqV2_eset }
\docType{data}
\title{ Integrated genomic analyses of ovarian carcinoma. }
\description{ A catalogue of molecular aberrations that cause ovarian cancer is critical for developing and deploying therapies that will improve patients' lives. The Cancer Genome Atlas project has analysed messenger RNA expression, microRNA expression, promoter methylation and DNA copy number in 489 high-grade serous ovarian adenocarcinomas and the DNA sequences of exons from coding genes in 316 of these tumours. Here we report that high-grade serous ovarian cancer is characterized by TP53 mutations in almost all tumours (96\%); low prevalence but statistically recurrent somatic mutations in nine further genes including NF1, BRCA1, BRCA2, RB1 and CDK12; 113 significant focal DNA copy number aberrations; and promoter methylation events involving 168 genes. Analyses delineated four ovarian cancer transcriptional subtypes, three microRNA subtypes, four promoter methylation subtypes and a transcriptional signature associated with survival duration, and shed new light on the impact that tumours with BRCA1/2 (BRCA1 or BRCA2) and CCNE1 aberrations have on survival. Pathway analyses suggested that homologous recombination is defective in about half of the tumours analysed, and that NOTCH and FOXM1 signalling are involved in serous ovarian cancer pathophysiology. }
\usage{data( TCGA.RNASeqV2_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Integrated genomic analyses of ovarian carcinoma. Nature 2011, 474:609-615. 
  Laboratory: Cancer Genome Atlas Research Network 2011 
  Contact information:  
  Title: Integrated genomic analyses of ovarian carcinoma. 
  URL:  
  PMIDs: 21720365 

  Abstract: A 179 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [RNASeqV2] Illumina HiSeq RNA sequencing
   platform_shorttitle:     
      Illumina HiSeq RNA sequencing
   platform_summary:     
      NA
   platform_manufacturer:     
      Illumina
   platform_distribution:     
      sequencing
   platform_accession:     
      NA
   platform_technology:     
      RNA sequencing

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ? A1BG ... ZZZ3 (20502 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 20502 features, 261 samples
Platform type: NA
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   5 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 256.00  256.00  256.00  143.00    3.62    3.19    4.03 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      261 character character 

unique_patient_ID: 
   Length     Class      Mode 
      261 character character 

sample_type: 
tumor 
  261 

histological_type: 
ser 
261 

primarysite: 
other    ov 
    1   260 

summarygrade: 
high  low NA's 
 226   29    6 

summarystage: 
early  late  NA's 
   18   242     1 

tumorstage: 
   2    3    4 NA's 
  18  209   33    1 

substage: 
   b    c NA's 
  16  211   34 

grade: 
   1    2    3    4 NA's 
   1   28  225    1    6 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  34.00   51.00   58.00   58.84   66.00   87.00 

pltx: 
   n    y NA's 
  17  215   29 

tax: 
   n    y NA's 
  17  215   29 

neo: 
   n NA's 
 232   29 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    9.0   225.0   426.5   585.3   755.0  5480.0      19 

recurrence_status: 
norecurrence   recurrence 
         123          138 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    9.0   341.8   878.0  1018.0  1446.0  5480.0       5 

vital_status: 
deceased   living     NA's 
     143      114        4 

site_of_tumor_first_recurrence: 
locoregional   metastasis         NA's 
          82           56          123 

primary_therapy_outcome_success: 
  completeresponse    partialresponse progressivedisease      stabledisease 
               147                 30                 15                 15 
              NA's 
                54 

debulking: 
   optimal suboptimal       NA's 
       171         60         30 

percent_normal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   0.000   0.000   2.066   0.000  55.000       5 

percent_stromal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    5.00   10.00   11.43   15.00   70.00       4 

percent_tumor_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00   77.00   85.00   82.07   90.00  100.00       4 

uncurated_author_metadata: 
   Length     Class      Mode 
      261 character character 

}}
\keyword{datasets}
