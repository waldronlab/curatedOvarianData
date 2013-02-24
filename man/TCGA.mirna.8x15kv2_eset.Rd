\name{ TCGA.mirna.8x15kv2_eset }
\alias{ TCGA.mirna.8x15kv2_eset }
\docType{data}
\title{ Integrated genomic analyses of ovarian carcinoma. }
\description{ A catalogue of molecular aberrations that cause ovarian cancer is critical for developing and deploying therapies that will improve patients' lives. The Cancer Genome Atlas project has analysed messenger RNA expression, microRNA expression, promoter methylation and DNA copy number in 489 high-grade serous ovarian adenocarcinomas and the DNA sequences of exons from coding genes in 316 of these tumours. Here we report that high-grade serous ovarian cancer is characterized by TP53 mutations in almost all tumours (96\%); low prevalence but statistically recurrent somatic mutations in nine further genes including NF1, BRCA1, BRCA2, RB1 and CDK12; 113 significant focal DNA copy number aberrations; and promoter methylation events involving 168 genes. Analyses delineated four ovarian cancer transcriptional subtypes, three microRNA subtypes, four promoter methylation subtypes and a transcriptional signature associated with survival duration, and shed new light on the impact that tumours with BRCA1/2 (BRCA1 or BRCA2) and CCNE1 aberrations have on survival. Pathway analyses suggested that homologous recombination is defective in about half of the tumours analysed, and that NOTCH and FOXM1 signalling are involved in serous ovarian cancer pathophysiology. }
\usage{data( TCGA.mirna.8x15kv2_eset )}
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
      [miRNA-8x15k2] Agilent Human miRNA G4470B
   platform_shorttitle:     
      Agilent miRNA-8x15k2 G4470B
   platform_summary:     
      NA
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      NA

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: ebv-miR-BART1-3p ebv-miR-BART1-5p ... kshv-miR-K12-9*
    (799 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 799 features, 554 samples
Platform type: NA
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   10 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 544.00  544.00  544.00  286.00    3.71    3.42    4.03 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      554 character character 

unique_patient_ID: 
   Length     Class      Mode 
      554 character character 

sample_type: 
tumor 
  554 

histological_type: 
ser 
554 

primarysite: 
other    ov 
    4   550 

summarygrade: 
high  low NA's 
 474   68   12 

summarystage: 
early  late  NA's 
   39   511     4 

tumorstage: 
   1    2    3    4 NA's 
  16   23  426   85    4 

substage: 
   b    c NA's 
  31  434   89 

grade: 
   1    2    3    4 NA's 
   6   62  473    1   12 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  26.00   51.00   59.00   59.81   69.00   89.00 

pltx: 
   n    y NA's 
  19  478   57 

tax: 
   n    y NA's 
  43  454   57 

neo: 
   n NA's 
 497   57 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   235.2   436.0   618.7   797.0  5480.0      44 

recurrence_status: 
norecurrence   recurrence 
         268          286 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   346.0   867.5   997.7  1446.0  5480.0      10 

vital_status: 
deceased   living     NA's 
     286      261        7 

site_of_tumor_first_recurrence: 
                locoregional locoregional_plus_metastatic 
                         145                            3 
                  metastasis                         NA's 
                         138                          268 

primary_therapy_outcome_success: 
  completeresponse    partialresponse progressivedisease      stabledisease 
               308                 63                 41                 30 
              NA's 
               112 

debulking: 
   optimal suboptimal       NA's 
       359        137         58 

percent_normal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   0.000   0.000   2.375   0.000  55.000      10 

percent_stromal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    5.00   10.00   12.78   19.00   70.00      16 

percent_tumor_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00   75.00   85.00   80.72   90.00  100.00      13 

uncurated_author_metadata: 
   Length     Class      Mode 
      554 character character 

}}
\keyword{datasets}
