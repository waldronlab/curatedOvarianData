\name{GSE6008_eset}
\alias{GSE6008_eset}
\docType{data}
\title{Human Ovarian Tumors and Normal Ovaries}
\description{Lysophosphatidic acid (LPA) governs a number of physiologic and pathophysiological processes. Malignant ascites fluid is rich in LPA, and LPA receptors are aberrantly expressed by ovarian cancer cells, implicating LPA in the initiation and progression of ovarian cancer. However, there is an absence of systematic data critically analyzing the transcriptional changes induced by LPA in ovarian cancer.}
\usage{data(GSE6008_eset)}
\format{

Authors: Murph MM, Liu W, Yu S, Lu Y, Hall H, Hennessy BT, Lahad J, Schaner M, Helland A, Kristensen G, Mills GB
  
Lab: Department of Systems Biology, University of Texas MD Anderson Cancer Center, Houston, Texas, United States of America.
  
Contact: mmurph@rx.uga.edu
  
Title: Lysophosphatidic acid-induced transcriptional profile represents serous epithelial ovarian carcinoma and worsened prognosis.

Url: http://pubget.com/search?q=19440550
  
PubMedID: 19440550
}
\details{
Platform used: [HG-U133A] Affymetrix Human Genome U133A Array.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 99 samples.

GEO_platform_accession: GSE6008

platform_summary: GPL96

biomart_ID: HG-U133A affy_hg_u133a

---------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------

alt_sample_name:

Length = 99.
Class =  character.
Mode = character.

sample_type:
tumor = 99.

subtype:
clearcell = 8. endo = 37. mucinous = 13. ser = 41.     

primarysite: 
ov = 99.

summarygrade:
high = 38.  low = 36. NA's = 25.

summarystage:
early = 42.  late = 53.  NA's = 4.

tumorstage
1 = 31.    2 = 11.    3 = 44.   4 = 9. NA's = 4.

substage:
a = 19.    b = 2.    c = 54.   d = 1. NA's = 23.

G:
1 = 19.    2 = 17.   3 = 38. NA's = 25.

debulking: 
unknown = 99.

batch:

2002-04-03 2002-04-04 2002-04-09 2002-04-10 2002-04-12 2002-08-13 2002-08-15 

         3          8          9          2          3          4          4 
         
2002-08-22 2002-08-23 2002-08-27 2002-08-28 2002-08-29 2002-08-30 2002-09-11 

         8          8          5          6         16         14          9 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE6008
}
\examples{
data(GSE6008_eset)
## maybe str(GSE6008_eset) ; plot(GSE6008_eset) ...
}
\keyword{datasets}
