\name{ GSE8842_eset }
\alias{ GSE8842_eset }
\docType{data}
\title{ Analysis of gene expression in early-stage ovarian cancer. }
\description{ Gene expression profile was analyzed in 68 stage I and 15 borderline ovarian cancers to determine if different clinical features of stage I ovarian cancer such as histotype, grade, and survival are related to differential gene expression.Tumors were obtained directly at surgery and immediately frozen in liquid nitrogen until analysis. Glass arrays containing 16,000 genes were used in a dual-color assay labeling protocol.Unsupervised analysis identified eight major patient partitions, one of which was statistically associated to overall survival, grading, and histotype and another with grading and histotype. Supervised analysis allowed detection of gene profiles clearly associated to histotype or to degree of differentiation. No difference was found between borderline and grade 1 tumors. As to recurrence, a subset of genes able to differentiate relapsers from nonrelapsers was identified. Among these, cyclin E and minichromosome maintenance protein 5 were found particularly relevant, as their expression was inversely correlated to progression-free survival (P = 0.00033 and 0.017, respectively).Specific molecular signatures define different histotypes and prognosis of stage I ovarian cancer. Mucinous and clear cells histotypes can be distinguished from the others regardless of tumor grade. Cyclin E and minichromosome maintenance protein 5, whose expression was found previously to be related to a bad prognosis of advanced ovarian cancer, appear to be potential prognostic markers in stage I ovarian cancer too, independent of other pathologic and clinical variables. }
\usage{data( GSE8842_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Marchini S, Mariani P, Chiorino G, Marrazzo E, Bonomi R, Fruscio R, Clivio L, Garbi A, Torri V, Cinquini M, Dell'Anna T, Apolone G, Broggini M, D'Incalci M.Analysis of gene expression in early-stage ovarian cancer. Clin Cancer Res. 2008 Dec 1; 14(23):7850-60. 
  Laboratory: Marchini, D'Incalci 2008 
  Contact information:  
  Title: Analysis of gene expression in early-stage ovarian cancer. 
  URL:  
  PMIDs: 19047114 

  Abstract: A 225 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent Human 1 cDNA Microarray (G4100A)
   platform_shorttitle:     
      Agilent  G4100A cDNA
   platform_summary:     
      hgug4100a
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      custom-commerical
   platform_accession:     
      GPL5689
   platform_technology:     
      spotted DNA/cDNA

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A2M AADAC ... ZYX (6536 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 6536 features, 83 samples
Platform type: hgug4100a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
     83      83      83      15      NA      12      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       83 character character 

sample_type: 
borderline      tumor 
        15         68 

histological_type: 
   Length     Class      Mode 
       83 character character 

primarysite: 
ov 
83 

summarygrade: 
high  low NA's 
  35   33   15 

summarystage: 
early 
   83 

tumorstage: 
 1 
83 

substage: 
 a  b  c 
25  5 53 

grade: 
   1    2    3 NA's 
  13   20   35   15 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  21.00   43.00   50.00   51.25   61.00   87.00 

recurrence_status: 
norecurrence   recurrence 
          62           21 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0    1192    2248    2273    3048    5824 

vital_status: 
deceased   living 
      15       68 

uncurated_author_metadata: 
   Length     Class      Mode 
       83 character character 

}}
\keyword{datasets}
