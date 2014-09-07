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
   platform_technology:     
      in situ oligonucleotide

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19816 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19816 features, 204 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      204 character character 

sample_type: 
    benign borderline      tumor 
         2          8        194 

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
   Length     Class      Mode 
      204 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      204 character character 

}}
\keyword{datasets}
