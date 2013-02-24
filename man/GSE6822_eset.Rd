\name{ GSE6822_eset }
\alias{ GSE6822_eset }
\docType{data}
\title{ Classification of ovarian tumor samples }
\description{ Ouellet V, Provencher DM, Maugard CM, Le Page C, Ren F, Lussier C, Novak J, Ge B, Hudson TJ, Tonin PN, Mes-Masson A-M: Discrimination between serous low malignant potential and invasive epithelial ovarian tumors using molecular profiling. Oncogene 2005, 24:4672-4687. }
\usage{data( GSE6822_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Ouellet V, Provencher DM, Maugard CM, Le Page C, Ren F, Lussier C, Novak J, Ge B, Hudson TJ, Tonin PN, Mes-Masson A-M: Discrimination between serous low malignant potential and invasive epithelial ovarian tumors using molecular profiling. Oncogene 2005, 24:4672-4687. 
  Laboratory: Ouellet, Mes-Masson 2005 
  Contact information:  
  Title: Classification of ovarian tumor samples 
  URL:  
  PMIDs: PMID unknown 

  Abstract: A 40 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [Hu6800] Affymetrix Human Full Length HuGeneFL Array
   platform_shorttitle:     
      Affymetrix Hu6800
   platform_summary:     
      hu6800
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL80

Preprocessing: rma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A2M AADAC ... ZYX (5231 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 5231 features, 66 samples
Platform type: hu6800
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       66 character character 

sample_type: 
tumor 
   66 

histological_type: 
   Length     Class      Mode 
       66 character character 

primarysite: 
ov 
66 

summarygrade: 
high  low NA's 
  40   15   11 

grade: 
   1    2    3 NA's 
   1   14   40   11 

batch: 
2000-12-21 2001-05-03 2001-05-29 2001-06-12 2001-09-25 2001-09-26 2001-09-27 
         1          1          3          3          1          5          8 
2002-02-14 2002-04-17 2002-04-18 2002-07-18 2002-07-24 2002-10-20 2002-10-30 
         4          1          9          7          4         10          5 
2002-11-01 2002-11-13 
         2          2 

uncurated_author_metadata: 
   Length     Class      Mode 
       66 character character 

}}
\keyword{datasets}
