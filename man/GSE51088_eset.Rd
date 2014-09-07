\name{ GSE51088_eset }
\alias{ GSE51088_eset }
\docType{data}
\title{ POSTN/TGFBI-associated stromal signature predicts poor prognosis in serous epithelial ovarian cancer. }
\description{ To identify molecular prognosticators and therapeutic targets for high-grade serous epithelial ovarian cancers (EOCs) using genetic analyses driven by biologic features of EOC pathogenesis.Ovarian tissue samples (n = 172; 122 serous EOCs, 30 other EOCs, 20 normal/benign) collected prospectively from sequential patients undergoing gynecologic surgery were analyzed using RNA expression microarrays. Samples were classified based on expression of genes with potential relevance in ovarian cancer. Gene sets were defined using Rosetta Similarity Search Tool (ROAST) and analysis of variance (ANOVA). Gene copy number variations were identified by array comparative genomic hybridization.No distinct subgroups of EOC could be identified by unsupervised clustering, however, analyses based on genes correlated with periostin (POSTN) and estrogen receptor-alpha (ESR1) yielded distinct subgroups. When 95 high-grade serous EOCs were grouped by genes based on ANOVA comparing ESR1/WT1 and POSTN/TGFBI samples, overall survival (OS) was significantly shorter for 43 patients with tumors expressing genes associated with POSTN/TGFBI compared to 52 patients with tumors expressing genes associated with ESR1/WT1 (median 30 versus 49 months, respectively; P = 0.022). Several targets with therapeutic potential were identified within each subgroup. BRCA germline mutations were more frequent in the ESR1/WT1 subgroup. Proliferation-associated genes and TP53 status (mutated or wild-type) did not correlate with survival. Findings were validated using independent ovarian cancer datasets.Two distinct molecular subgroups of high-grade serous EOCs based on POSTN/TGFBI and ESR1/WT1 expressions were identified with significantly different OS. Specific differentially expressed genes between these subgroups provide potential prognostic and therapeutic targets.Copyright ?? 2013 Elsevier Inc. All rights reserved. }
\usage{data( GSE51088_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Karlan BY, Dering J, Walsh C, Orsulic S, Lester J, Anderson LA, Ginther CL, Fejzo M, Slamon D 
  Laboratory: Karlan, Slamon 2014 
  Contact information:  
  Title: POSTN/TGFBI-associated stromal signature predicts poor prognosis in serous epithelial ovarian cancer. 
  URL:  
  PMIDs: 24368280 

  Abstract: A 250 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent-012097 Human 1A Microarray (V2) G4110B (Probe Name version)
   platform_shorttitle:     
      Agilent G4110B
   platform_summary:     
      hgug4110b
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      GPL7264
   platform_technology:     
      in situ oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (8211 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 8211 features, 172 samples
Platform type: hgug4110b
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   20 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 152.00  152.00  152.00  112.00    4.13    3.50    4.92 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      172 character character 

sample_type: 
    benign borderline    healthy metastatic      tumor 
         5         12         15         17        123 

histological_type: 
   Length     Class      Mode 
      172 character character 

summarygrade: 
high  low NA's 
 119   30   23 

summarystage: 
early  late  NA's 
   31   120    21 

tumorstage: 
   1    2    3    4 NA's 
  22    9  103   17   21 

substage: 
   a    b    c NA's 
  17   22   94   39 

grade: 
   0    1    2    3 NA's 
   8    8   14  119   23 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   26.0    49.0    57.5    58.6    68.0    91.0 

neo: 
  n 
172 

recurrence_status: 
norecurrence   recurrence         NA's 
          36          111           25 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
     30     791    1491    1835    2344    7001      20 

vital_status: 
deceased   living     NA's 
     112       40       20 

percent_normal_cells: 
 30- NA's 
 140   32 

percent_stromal_cells: 
 30- NA's 
 140   32 

percent_tumor_cells: 
 70+ NA's 
 140   32 

uncurated_author_metadata: 
   Length     Class      Mode 
      172 character character 

}}
\keyword{datasets}
