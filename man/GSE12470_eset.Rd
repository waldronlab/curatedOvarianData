\name{ GSE12470_eset }
\alias{ GSE12470_eset }
\docType{data}
\title{ Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis. }
\description{ To elucidate the mechanisms of rapid progression of serous ovarian cancer, gene expression profiles from 43 ovarian cancer tissues comprising eight early stage and 35 advanced stage tissues were carried out using oligonucleotide microarrays of 18,716 genes. By non-negative matrix factorization analysis using 178 genes, which were extracted as stage-specific genes, 35 advanced stage cases were classified into two subclasses with superior (n = 17) and poor (n = 18) outcome evaluated by progression-free survival (log rank test, P = 0.03). Of the 178 stage-specific genes, 112 genes were identified as showing different expression between the two subclasses. Of the 48 genes selected for biological function by gene ontology analysis or Ingenuity Pathway Analysis, five genes (ZEB2, CDH1, LTBP2, COL16A1, and ACTA2) were extracted as candidates for prognostic factors associated with progression-free survival. The relationship between high ZEB2 or low CDH1 expression and shorter progression-free survival was validated by real-time RT-PCR experiments of 37 independent advanced stage cancer samples. ZEB2 expression was negatively correlated with CDH1 expression in advanced stage samples, whereas ZEB2 knockdown in ovarian adenocarcinoma SKOV3 cells resulted in an increase in CDH1 expression. Multivariate analysis showed that high ZEB2 expression was independently associated with poor prognosis. Furthermore, the prognostic effect of E-cadherin encoded by CDH1 was verified using immunohistochemical analysis of an independent advanced stage cancer samples set (n = 74). These findings suggest that the expression of epithelial-mesenchymal transition-related genes such as ZEB2 and CDH1 may play important roles in the invasion process of advanced stage serous ovarian cancer. }
\usage{data( GSE12470_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Yoshihara K, Tajima A, Komata D, Yamamoto T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo Y, Inoue I, Tanaka K.Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis. Cancer Sci. 2009 Aug; 100(8):1421-8. 
  Laboratory: Yoshihara, Tanaka 2009 
  Contact information:  
  Title: Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis. 
  URL:  
  PMIDs: 19486012 

  Abstract: A 253 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      Agilent-012097 Human 1A Microarray (V2) G4110B (Feature Number version)
   platform_shorttitle:     
      Agilent G4110B
   platform_summary:     
      hgug4110b
   platform_manufacturer:     
      Agilent
   platform_distribution:     
      commercial
   platform_accession:     
      GPL887
   platform_technology:     
      in situ oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (16889 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 16889 features, 53 samples
Platform type: hgug4110b
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       53 character character 

sample_type: 
healthy   tumor 
     10      43 

histological_type: 
 ser NA's 
  43   10 

primarysite: 
ov 
53 

summarystage: 
early  late  NA's 
    8    35    10 

tumorstage: 
   1 NA's 
   8   45 

uncurated_author_metadata: 
   Length     Class      Mode 
       53 character character 

}}
\keyword{datasets}
