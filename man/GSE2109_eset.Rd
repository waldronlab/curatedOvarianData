\name{ GSE2109_eset }
\alias{ GSE2109_eset }
\docType{data}
\title{ IGC EXpression Project for Oncology }
\description{ EXpression Project for Oncology, International Genomics Consortium, www.intgen.org }
\usage{data( GSE2109_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: EXpression Project for Oncology, International Genomics Consortium, www.intgen.org 
  Laboratory: expO, IGC 2005 
  Contact information:  
  Title: IGC EXpression Project for Oncology 
  URL:  
  PMIDs: PMID unknown 

  Abstract: A 8 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      Affymetrix HG-U133Plus2
   platform_summary:     
      hgu133plus2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL570

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (19093 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19093 features, 204 samples
Platform type: hgu133plus2
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
   Length     Class      Mode 
      204 character character 

primarysite: 
other    ov  NA's 
   23   178     3 

summarygrade: 
high  low NA's 
  91   31   82 

summarystage: 
early  late  NA's 
   37    87    80 

tumorstage: 
   1    2    3    4 NA's 
  20   14   58   18   94 

substage: 
   a    b    c NA's 
  17   22   79   86 

grade: 
   1    2    3    4 NA's 
  11   20   83    8   82 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  25.00   45.00   55.00   58.82   65.00   85.00 

batch: 
2004-12-03 2004-12-04 2004-12-07 2005-02-11 2005-03-03 2005-03-10 2005-03-11 
         3          3          1          1          1          1          1 
2005-03-15 2005-03-16 2005-03-17 2005-03-19 2005-03-22 2005-04-13 2005-04-26 
         3          1          2          4          2          1          5 
2005-04-29 2005-05-10 2005-06-01 2005-06-03 2005-06-08 2005-06-17 2005-08-05 
         2          2          5          3          3          6          3 
2005-08-09 2005-08-11 2005-09-07 2005-09-09 2005-09-13 2005-11-02 2005-11-04 
         1          6          1          3          3          6          3 
2005-11-15 2005-11-18 2005-12-02 2006-01-24 2006-01-26 2006-02-07 2006-02-28 
         3          1          4          2          1          1          1 
2006-03-06 2006-03-14 2006-04-18 2006-04-20 2006-05-16 2006-06-08 2006-07-26 
         2          2          1          2          3          1          2 
2006-07-28 2006-09-12 2006-09-14 2006-10-10 2006-10-24 2006-10-31 2006-11-09 
         1          2          1          1          9          5         10 
2006-11-21 2006-11-30 2006-12-07 2007-01-12 2007-02-09 2007-03-07 2007-03-09 
         1          6          3          1          1          8          1 
2007-03-15 2007-05-01 2007-05-03 2007-05-15 2007-05-18 2007-05-30 2007-06-12 
         4          2          3          4          2          2          1 
2007-07-27 2007-09-05 2007-09-07 2007-09-11 2007-09-12 2008-02-15 2008-02-21 
         2          3          1          4          4          1          3 
2008-02-27 2008-03-04 2008-05-13 2008-05-16 2008-05-23 
         2          1          4          4          5 

uncurated_author_metadata: 
   Length     Class      Mode 
      204 character character 

}}
\keyword{datasets}
