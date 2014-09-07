\name{ E.MTAB.386_eset }
\alias{ E.MTAB.386_eset }
\docType{data}
\title{ Angiogenic mRNA and microRNA gene expression signature predicts a novel subtype of serous ovarian cancer. }
\description{ Ovarian cancer is the fifth leading cause of cancer death for women in the U.S. and the seventh most fatal worldwide. Although ovarian cancer is notable for its initial sensitivity to platinum-based therapies, the vast majority of patients eventually develop recurrent cancer and succumb to increasingly platinum-resistant disease. Modern, targeted cancer drugs intervene in cell signaling, and identifying key disease mechanisms and pathways would greatly advance our treatment abilities. In order to shed light on the molecular diversity of ovarian cancer, we performed comprehensive transcriptional profiling on 129 advanced stage, high grade serous ovarian cancers. We implemented a, re-sampling based version of the ISIS class discovery algorithm (rISIS: robust ISIS) and applied it to the entire set of ovarian cancer transcriptional profiles. rISIS identified a previously undescribed patient stratification, further supported by micro-RNA expression profiles, and gene set enrichment analysis found strong biological support for the stratification by extracellular matrix, cell adhesion, and angiogenesis genes. The corresponding "angiogenesis signature" was validated in ten published independent ovarian cancer gene expression datasets and is significantly associated with overall survival. The subtypes we have defined are of potential translational interest as they may be relevant for identifying patients who may benefit from the addition of anti-angiogenic therapies that are now being tested in clinical trials. }
\usage{data( E.MTAB.386_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Bentink S, Haibe-Kains B, Risch T, Fan J-B, Hirsch MS, Holton K, Rubio R, April C, Chen J, Wickham-Garcia E, Liu J, Culhane A, Drapkin R, Quackenbush J, Matulonis UA: Angiogenic mRNA and microRNA gene expression signature predicts a novel subtype of serous ovarian cancer. PLoS ONE 2012, 7:e30269. 
  Laboratory: Bentink, Matulonis 2012 
  Contact information:  
  Title: Angiogenic mRNA and microRNA gene expression signature predicts a novel subtype of serous ovarian cancer. 
  URL:  
  PMIDs: 22348002 

  Abstract: A 212 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Illumina humanRef-8 v2.0 expression beadchip
   platform_shorttitle:     
      Illumina humanRef-8 v2.0
   platform_summary:     
      illuminaHumanv2
   platform_manufacturer:     
      Illumina
   platform_distribution:     
      commercial
   platform_accession:     
      GPL6104
   platform_technology:     
      in situ oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A2M A4GALT ... ZZEF1 (10357 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 10357 features, 129 samples
Platform type: illuminaHumanv2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 129.00  129.00  129.00   73.00    3.51    2.68    4.13 

--------------------------- 
Available sample meta-data: 
--------------------------- 

unique_patient_ID: 
   Length     Class      Mode 
      129 character character 

sample_type: 
tumor 
  129 

histological_type: 
ser 
129 

primarysite: 
 ov 
129 

summarygrade: 
high 
 129 

summarystage: 
early  late 
    1   128 

tumorstage: 
  2   3   4 
  1 109  19 

substage: 
   a    b    c NA's 
   5   12   93   19 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   50.00   66.00   60.71   72.00   95.00 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    3.9   516.9   917.1  1007.0  1401.0  2724.0 

vital_status: 
deceased   living 
      73       56 

debulking: 
   optimal suboptimal       NA's 
        98         28          3 

uncurated_author_metadata: 
   Length     Class      Mode 
      129 character character 

}}
\keyword{datasets}
