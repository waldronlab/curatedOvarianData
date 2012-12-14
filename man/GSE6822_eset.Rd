\name{GSE6822_eset}
\alias{GSE6822_eset}
\docType{data}
\title{Classification of ovarian tumor samples}
\description{We have previously described differences in the expression profiles of tumors by comparing tumors of two known classes: low malignant potential versus highly malignant (invasive). This report presents an effort to discover unknown classes within a large heterogeneous set of ovarian tumors, using unsupervised learning approaches.  We were able to define four classes in a set of 74 ovarian cancer tumors. Groups identified as A1 and A2 were closely related and correlated with the invasive pathology class while the group identified as B2 was correlated with low malignant potential tumors; group B1 consisted of a mixture of low potential and invasive samples. We selected characteristic candidate genes, which were validated by quantitative PCR and by comparison to other published studies.}
\usage{data(GSE6822_eset)}
\format{

Authors: Novak JP, Ouellet V, Le Page C, Martinu K, Ponton A, Bachvarow DD, Filali-Mouhim A, Provencher DM, Tonin PN, Hudson TJ, Mes-Masson A.

Lab: McGill University and Genome Quebec Innovation Centre.

Contact: jaroslav.novak@mail.mcgill.ca, jaroslav.novak@gmail.com

Title: Classification of ovarian tumor samples.

Url: http://www.genomequebec.mcgill.ca/ovarian
}
\details{
Platform used: [Hu6800] Affymetrix Human Full Length HuGeneFL Array.

rma normalization used from the rma bioconductor package.

assayData: 5528 features, 66 samples.

GEO_platform_accession: GSE6822

platform_summary: GPL80

biomart_ID: Hu6800 affy_hugenefl

----------------------------------------------------------------------


Available sample meta-data: 

----------------------------------------------------------------------

alt_sample_name:
Length = 66.
Class =  character.
Mode = character.

sample_type: 
tumor = 66.

subtype:
Length = 66.
Class = character.
Mode = character.

primarysite: 
ov = 66.

summarygrade:
high = 40. low = 15. NA's = 11.

G:
1 = 1.  2 = 14.    3 = 40. NA's = 11.

debulking: 
unknown = 66.

batch: 
2000-12-21 2001-05-03 2001-05-29 2001-06-12 2001-09-25 2001-09-26 2001-09-27 

         1          1          3          3          1          5          8 
         
2002-02-14 2002-04-17 2002-04-18 2002-07-18 2002-07-24 2002-10-20 2002-10-30 

         4          1          9          7          4         10          5 
         
2002-11-01 2002-11-13 
         2          2 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE6822
}
\examples{
data(GSE6822_eset)
## maybe str(GSE6822_eset) ; plot(GSE6822_eset) ...
}
\keyword{datasets}
