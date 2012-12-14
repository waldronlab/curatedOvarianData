\name{GSE20565_eset}
\alias{GSE20565_eset}
\docType{data}
\title{A genomic and transcriptomic approach for a differential diagnosis between primary and secondary ovarian carcinomas in patients with a previous history of breast cancer.}
\description{The distinction between primary and secondary ovarian tumors may be challenging for pathologists. The purpose of the present work was to develop genomic and transcriptomic tools to further refine the pathological diagnosis of ovarian tumors after a previous history of breast cancer.}
\usage{data(GSE20565_eset)}
\format{

Authors: Meyniel JP, Cottu PH, Decraene C, Stern MH, Couturier J, Lebigot I, Nicolas A, Weber N, Fourchotte V, Alran S, Rapinat A, Gentien D, Roman-Roman S, Mignot L, Sastre-Garau X.

Lab: Department of Translational Research, Institut Curie, 26 rue d'Ulm, 75248 Paris, Cedex 05, France.

Contact: jean-philippe.meyniel@curie.fr

Title: A genomic and transcriptomic approach for a differential diagnosis between primary and secondary ovarian carcinomas in patients with a previous history of breast cancer.

Url: http://pubget.com/search?q=20492709.

PubMedID: 20492709.
 
}
\details{
Platforms used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0.

In the study but not included in curatedOvarianData:

   **Array|[Mapping50K_Xba240] Affymetrix Human Mapping 50K Xba240 SNP.

   **Array|[GenomeWideSNP_6] Affymetrix Genome-Wide Human SNP 6.0 Array.

frma normalization used from frma bioconductor package.

assayData: 18708 features, 140 samples.

-------------------------------------------------------------------------- 


Available sample meta-data: 

--------------------------------------------------------------------------

alt_sample_name: 
Length = 140.
Class = character.
Mode = character.

sample_type: 
tumor = 140.

subtype: 
Length = 140.
Class = character.
Mode = character.

primarysite: 
other = 44. ov = 96.
   
summarygrade: 
high = 63. low = 33. NA's = 44.

summarystage: 
early = 27. late = 67. NA's = 46.

tumorstage:
1 = 18. 2 = 9. 3 = 52. 4  = 15. NA's = 46.

substage:
a = 14. b = 10. c = 55. NA's = 61. 

G:
1 = 6. 2  = 27. 3 = 63. NA's = 44.

M:
0 = 96.  1 = 44. 

debulking:
unknown = 140.

batch: 

2006-06-01 2006-06-27 2006-06-28 2006-06-29 2006-06-30 2006-07-20 2008-03-06 

        21         18         37         20         36          7          1 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE20565
}
\examples{
data(GSE20565_eset)
## maybe str(GSE20565_eset) ; plot(GSE20565_eset) ...
}
\keyword{datasets}
