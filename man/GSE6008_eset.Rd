\name{ GSE6008_eset }
\alias{ GSE6008_eset }
\docType{data}
\title{ Lysophosphatidic acid-induced transcriptional profile represents serous epithelial ovarian carcinoma and worsened prognosis. }
\description{ Lysophosphatidic acid (LPA) governs a number of physiologic and pathophysiological processes. Malignant ascites fluid is rich in LPA, and LPA receptors are aberrantly expressed by ovarian cancer cells, implicating LPA in the initiation and progression of ovarian cancer. However, there is an absence of systematic data critically analyzing the transcriptional changes induced by LPA in ovarian cancer.In this study, gene expression profiling was used to examine LPA-mediated transcription by exogenously adding LPA to human epithelial ovarian cancer cells for 24 h to mimic long-term stimulation in the tumor microenvironment. The resultant transcriptional profile comprised a 39-gene signature that closely correlated to serous epithelial ovarian carcinoma. Hierarchical clustering of ovarian cancer patient specimens demonstrated that the signature is associated with worsened prognosis. Patients with LPA-signature-positive ovarian tumors have reduced disease-specific and progression-free survival times. They have a higher frequency of stage IIIc serous carcinoma and a greater proportion is deceased. Among the 39-gene signature, a group of seven genes associated with cell adhesion recapitulated the results. Out of those seven, claudin-1, an adhesion molecule and phenotypic epithelial marker, is the only independent biomarker of serous epithelial ovarian carcinoma. Knockdown of claudin-1 expression in ovarian cancer cells reduces LPA-mediated cellular adhesion, enhances suspended cells and reduces LPA-mediated migration.The data suggest that transcriptional events mediated by LPA in the tumor microenvironment influence tumor progression through modulation of cell adhesion molecules like claudin-1 and, for the first time, report an LPA-mediated expression signature in ovarian cancer that predicts a worse prognosis. }
\usage{data( GSE6008_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Murph MM, Liu W, Yu S, Lu Y, Hall H, Hennessy BT, Lahad J, Schaner M, Helland A, Kristensen G, Brresen-Dale AL, Mills GB. Lysophosphatidic acid-induced transcriptional profile represents serous epithelial ovarian carcinoma and worsened prognosis. PLoS One. 2009; 4(5):e5583. 
  Laboratory: Murph, Mills 2009 
  Contact information:  
  Title: Lysophosphatidic acid-induced transcriptional profile represents serous epithelial ovarian carcinoma and worsened prognosis. 
  URL:  
  PMIDs: 19440550 

  Abstract: A 247 word abstract is available. Use 'abstract' method.
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

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12981 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12981 features, 103 samples
Platform type: hgu133a
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      103 character character 

sample_type: 
healthy   tumor 
      4      99 

histological_type: 
clearcell      endo  mucinous       ser      NA's 
        8        37        13        41         4 

primarysite: 
 ov 
103 

summarygrade: 
high  low NA's 
  38   36   29 

summarystage: 
early  late  NA's 
   42    53     8 

tumorstage: 
   1    2    3    4 NA's 
  35   11   44    9    4 

substage: 
   a    b    c    d NA's 
  19    2   54    1   27 

grade: 
   1    2    3 NA's 
  19   17   38   29 

batch: 
2002-04-03 2002-04-04 2002-04-09 2002-04-10 2002-04-12 2002-08-13 2002-08-15 
         3          8          9          2          3          4          4 
2002-08-22 2002-08-23 2002-08-27 2002-08-28 2002-08-29 2002-08-30 2002-09-11 
         8          8          5          6         16         14          9 
2006-01-27 
         4 

uncurated_author_metadata: 
   Length     Class      Mode 
      103 character character 

}}
\keyword{datasets}
