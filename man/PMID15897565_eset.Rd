\name{ PMID15897565_eset }
\alias{ PMID15897565_eset }
\docType{data}
\title{ Patterns of gene expression that characterize long-term survival in advanced stage serous ovarian cancers. }
\description{ A better understanding of the underlying biology of invasive serous ovarian cancer is critical for the development of early detection strategies and new therapeutics. The objective of this study was to define gene expression patterns associated with favorable survival.RNA from 65 serous ovarian cancers was analyzed using Affymetrix U133A microarrays. This included 54 stage III/IV cases (30 short-term survivors who lived <3 years and 24 long-term survivors who lived >7 years) and 11 stage I/II cases. Genes were screened on the basis of their level of and variability in expression, leaving 7,821 for use in developing a predictive model for survival. A composite predictive model was developed that combines Bayesian classification tree and multivariate discriminant models. Leave-one-out cross-validation was used to select and evaluate models.Patterns of genes were identified that distinguish short-term and long-term ovarian cancer survivors. The expression model developed for advanced stage disease classified all 11 early-stage ovarian cancers as long-term survivors. The MAL gene, which has been shown to confer resistance to cancer therapy, was most highly overexpressed in short-term survivors (3-fold compared with long-term survivors, and 29-fold compared with early-stage cases). These results suggest that gene expression patterns underlie differences in outcome, and an examination of the genes that provide this discrimination reveals that many are implicated in processes that define the malignant phenotype.Differences in survival of advanced ovarian cancers are reflected by distinct patterns of gene expression. This biological distinction is further emphasized by the finding that early-stage cancers share expression patterns with the advanced stage long-term survivors, suggesting a shared favorable biology. }
\usage{data( PMID15897565_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Berchuck A, Iversen ES, Lancaster JM, Pittman J, Luo J, Lee P, Murphy S, Dressman HK, Febbo PG, West M, Nevins JR, Marks JR.Patterns of gene expression that characterize long-term survival in advanced stage serous ovarian cancers. Clin Cancer Res. 2005 May 15; 11(10):3686-96. 
  Laboratory: Berchuck, Marks 2005 
  Contact information:  
  Title: Patterns of gene expression that characterize long-term survival in advanced stage serous ovarian cancers. 
  URL:  
  PMIDs: 15897565 

  Abstract: A 258 word abstract is available. Use 'abstract' method.
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
      These samples are a subset of PMID17290060.

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12981 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12981 features, 63 samples
Platform type: hgu133a
Binary overall survival summary (definitions of long and short provided by study authors): 

 long short  NA's 
   24    28    11 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   1761    1828    1907    2001    2032    2536 

sample_type: 
tumor 
   63 

histological_type: 
ser 
 63 

primarysite: 
ov 
63 

summarygrade: 
high  low NA's 
  25   37    1 

summarystage: 
early  late 
   11    52 

tumorstage: 
 1  2  3  4 
 7  4 48  4 

grade: 
   1    2    3    4 NA's 
   2   35   24    1    1 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  33.00   52.50   59.00   59.21   67.00   79.00 

os_binary: 
 long short  NA's 
   24    28    11 

debulking: 
   optimal suboptimal       NA's 
        24         28         11 

batch: 
2002-09-20 2002-10-23 2002-11-12 2002-12-16 2002-12-21 2003-01-03 2003-05-30 
        15          9         10          1          3         11         13 
2003-07-02 
         1 

uncurated_author_metadata: 
   Length     Class      Mode 
       63 character character 

}}
\keyword{datasets}
