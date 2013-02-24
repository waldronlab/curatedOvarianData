\name{ PMID17290060_eset }
\alias{ PMID17290060_eset }
\docType{data}
\title{ An integrated genomic-based approach to individualized treatment of patients with advanced-stage ovarian cancer. }
\description{ The purpose of this study was to develop an integrated genomic-based approach to personalized treatment of patients with advanced-stage ovarian cancer. We have used gene expression profiles to identify patients likely to be resistant to primary platinum-based chemotherapy and also to identify alternate targeted therapeutic options for patients with de novo platinum-resistant disease.A gene expression model that predicts response to platinum-based therapy was developed using a training set of 83 advanced-stage serous ovarian cancers and tested on a 36-sample external validation set. In parallel, expression signatures that define the status of oncogenic signaling pathways were evaluated in 119 primary ovarian cancers and 12 ovarian cancer cell lines. In an effort to increase chemotherapy sensitivity, pathways shown to be activated in platinum-resistant cancers were subject to targeted therapy in ovarian cancer cell lines.Gene expression profiles identified patients with ovarian cancer likely to be resistant to primary platinum-based chemotherapy with greater than 80\% accuracy. In patients with platinum-resistant disease, we identified expression signatures consistent with activation of Src and Rb/E2F pathways, components of which were successfully targeted to increase response in ovarian cancer cell lines.We have defined a strategy for treatment of patients with advanced-stage ovarian cancer that uses therapeutic stratification based on predictions of response to chemotherapy, coupled with prediction of oncogenic pathway deregulation, as a method to direct the use of targeted agents. }
\usage{data( PMID17290060_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Dressman HK, Berchuck A, Chan G, Zhai J, Bild A, Sayer R, Cragun J, Clarke J, Whitaker RS, Li L, Gray J, Marks J, Ginsburg GS, Potti A, West M, Nevins JR, Lancaster JM.An integrated genomic-based approach to individualized treatment of patients with advanced-stage ovarian cancer. J Clin Oncol. 2007 Feb 10; 25(5):517-25. 
  Laboratory: Dressman, Lancaster 2007 
  Contact information:  
  Title: An integrated genomic-based approach to individualized treatment of patients with advanced-stage ovarian cancer. 
  URL:  
  PMIDs: 17290060 

  Abstract: A 223 word abstract is available. Use 'abstract' method.
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
      This paper has been retracted.

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12981 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12981 features, 117 samples
Platform type: hgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 117.00  117.00  117.00   67.00    5.26    2.79    7.48 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      117 character character 

sample_type: 
tumor 
  117 

histological_type: 
ser 
117 

primarysite: 
 ov 
117 

summarygrade: 
high  low NA's 
  57   57    3 

summarystage: 
early  late  NA's 
    1   115     1 

tumorstage: 
   2    3    4 NA's 
   1   98   17    1 

grade: 
   1    2    3    4 NA's 
   4   53   56    1    3 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
     30     510    1020    1496    2220    5550 

vital_status: 
deceased   living 
      67       50 

primary_therapy_outcome_success: 
  completeresponse progressivedisease 
                85                 32 

debulking: 
   optimal suboptimal 
        63         54 

batch: 
2002-09-20 2002-10-23 2002-11-12 2002-12-16 2002-12-21 2003-01-03 2003-05-30 
        10          8          9          1          3         11         10 
2004-03-09 2004-03-16 2004-04-20 2004-05-18 2004-05-21 2004-05-27 2004-06-22 
        16          6          5         15          7          7          1 
2004-06-23 
         8 

uncurated_author_metadata: 
   Length     Class      Mode 
      117 character character 

}}
\keyword{datasets}
