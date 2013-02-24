\name{ TCGA_eset }
\alias{ TCGA_eset }
\docType{data}
\title{ Integrated genomic analyses of ovarian carcinoma. }
\description{ A catalogue of molecular aberrations that cause ovarian cancer is critical for developing and deploying therapies that will improve patients' lives. The Cancer Genome Atlas project has analysed messenger RNA expression, microRNA expression, promoter methylation and DNA copy number in 489 high-grade serous ovarian adenocarcinomas and the DNA sequences of exons from coding genes in 316 of these tumours. Here we report that high-grade serous ovarian cancer is characterized by TP53 mutations in almost all tumours (96\%); low prevalence but statistically recurrent somatic mutations in nine further genes including NF1, BRCA1, BRCA2, RB1 and CDK12; 113 significant focal DNA copy number aberrations; and promoter methylation events involving 168 genes. Analyses delineated four ovarian cancer transcriptional subtypes, three microRNA subtypes, four promoter methylation subtypes and a transcriptional signature associated with survival duration, and shed new light on the impact that tumours with BRCA1/2 (BRCA1 or BRCA2) and CCNE1 aberrations have on survival. Pathway analyses suggested that homologous recombination is defective in about half of the tumours analysed, and that NOTCH and FOXM1 signalling are involved in serous ovarian cancer pathophysiology. }
\usage{data( TCGA_eset )}
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
      [HT_HG-U133A] Affymetrix HT Human Genome U133A Array
   platform_shorttitle:     
      Affymetrix HT_HG-U133A
   platform_summary:     
      hthgu133a
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL3921
   warnings:     
      The following samples are likely from specimens also used in GSE26712: TCG
A.13.0725, TCGA.13.0885, TCGA.13.0887, TCGA.13.0890, TCGA.13.0886, TCGA.13
.0714, TCGA.13.0727, TCGA.13.1817, TCGA.13.1499, TCGA.13.0883

Preprocessing: rma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12981 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12981 features, 578 samples
Platform type: hthgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   21 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 557.00  557.00  557.00  290.00    3.73    3.45    4.06 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      578 character character 

unique_patient_ID: 
   Length     Class      Mode 
      578 character character 

sample_type: 
adjacentnormal          tumor 
             8            570 

histological_type: 
 ser NA's 
 568   10 

primarysite: 
other    ov  NA's 
    4   564    10 

summarygrade: 
high  low NA's 
 480   75   23 

summarystage: 
early  late  NA's 
   43   520    15 

tumorstage: 
   1    2    3    4 NA's 
  16   27  436   84   15 

substage: 
   b    c NA's 
  31  448   99 

grade: 
   1    2    3    4 NA's 
   6   69  479    1   23 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   51.00   59.00   59.70   68.25   89.00      10 

pltx: 
   n    y NA's 
  19  492   67 

tax: 
   n    y NA's 
  43  468   67 

neo: 
   n NA's 
 511   67 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    8.0   238.2   443.5   623.7   812.0  5480.0      56 

recurrence_status: 
norecurrence   recurrence 
         279          299 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
      8     349     881    1010    1446    5480      21 

vital_status: 
deceased   living     NA's 
     290      270       18 

site_of_tumor_first_recurrence: 
                locoregional locoregional_plus_metastatic 
                         153                            3 
                  metastasis                         NA's 
                         143                          279 

primary_therapy_outcome_success: 
  completeresponse    partialresponse progressivedisease      stabledisease 
               318                 65                 41                 30 
              NA's 
               124 

debulking: 
   optimal suboptimal       NA's 
       367        140         71 

percent_normal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.000   0.000   0.000   2.385   0.000  55.000      19 

percent_stromal_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    5.00   10.00   12.85   20.00   70.00      25 

percent_tumor_cells: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00   75.00   85.00   80.64   90.00  100.00      22 

batch: 
   9   11   12   13   14   15   17   18   19   21   22   24   27   40 NA's 
  45   37   45   47   46   22   47   47   47   46   47   46    6   49    1 

uncurated_author_metadata: 
   Length     Class      Mode 
      578 character character 

}}
\keyword{datasets}
