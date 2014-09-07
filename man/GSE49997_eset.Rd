\name{ GSE49997_eset }
\alias{ GSE49997_eset }
\docType{data}
\title{ Validating the impact of a molecular subtype in ovarian cancer on outcomes: a study of the OVCAD Consortium. }
\description{ Most patients with epithelial ovarian cancer (EOC) are diagnosed at advanced stage and have a poor prognosis. However, a small proportion of these patients will survive, whereas others will die very quickly. Clinicopathological factors do not allow precise identification of these subgroups. Thus, we have validated a molecular subclassification as new prognostic factor in EOC. One hundred and ninety-four patients with Stage II-IV EOC were characterized by whole-genome expression profiling of tumor tissues and were classified using a published 112 gene set, derived from an International Federation of Gynecology and Obstetrics (FIGO) stage-directed supervised classification approach. The 194 tumor samples were classified into two subclasses comprising 95 (Subclass 1) and 99 (Subclass 2) tumors. All nine FIGO II tumors were grouped in Subclass 1 (P = 0.001). Subclass 2 (54\% of advanced-stage tumors) was significantly correlated with peritoneal carcinomatosis and non-optimal debulking. Patients with Subclass 2 tumors had a worse overall survival for both serous and non-serous histological subtypes, as revealed by univariate analysis (hazard ratios [HR] of 3.17 and 17.11, respectively; P ??? 0.001) and in models corrected for relevant clinicopathologic parameters (HR 2.87 and 12.42, respectively; P ??? 0.023). Significance analysis of microarrays revealed 2082 genes that were differentially expressed in advanced-grade serous tumors of both subclasses and the focal adhesion pathway as the most deregulated pathway. In the present validation study, we have shown that, in advanced-stage serous ovarian cancer, two approximately equally large molecular subtypes exist, independent of classical clinocopathological parameters and presenting with highly different whole-genome expression profiles and a markedly different overall survival. Similar results were obtained in a small cohort of patients with non-serous tumors.?? 2012 Japanese Cancer Association. }
\usage{data( GSE49997_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Pils D1, Hager G, Tong D, Aust S, Heinze G, Kohl M, Schuster E, Wolf A, Sehouli J, Braicu I, Vergote I, Cadron I, Mahner S, Hofstetter G, Speiser P, Zeillinger R 
  Laboratory: Pils, Zeilinger 2012 
  Contact information:  
  Title: Validating the impact of a molecular subtype in ovarian cancer on outcomes: a study of the OVCAD Consortium. 
  URL:  
  PMIDs: 22497737 

  Abstract: A 276 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      ABI Human Genome Survey Microarray Version 2
   platform_shorttitle:     
      ABI Human Genome
   platform_summary:     
      
   platform_manufacturer:     
      Applied Biosystems
   platform_distribution:     
      commercial
   platform_accession:     
      GPL2986
   platform_technology:     
      in situ oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (16048 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 16048 features, 204 samples
Platform type: 
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   10 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 194.00  194.00  194.00   57.00      NA    3.67      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      204 character character 

sample_type: 
tumor 
  204 

histological_type: 
other   ser  NA's 
   23   171    10 

summarygrade: 
high  low NA's 
 143   50   11 

summarystage: 
early  late  NA's 
    9   185    10 

tumorstage: 
   2    3    4 NA's 
   9  154   31   10 

grade: 
   2    3 NA's 
  50  143   11 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   50.00   57.00   57.66   67.00   85.00      10 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   30.0   335.0   487.0   580.1   722.5  1461.0      10 

recurrence_status: 
norecurrence   recurrence         NA's 
          70          124           10 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   30.0   517.0   745.5   782.9  1027.0  1491.0      10 

vital_status: 
deceased   living     NA's 
      57      137       10 

debulking: 
   optimal suboptimal       NA's 
       137         57         10 

uncurated_author_metadata: 
   Length     Class      Mode 
      204 character character 

}}
\keyword{datasets}
