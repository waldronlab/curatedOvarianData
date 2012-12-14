\name{GSE12470_eset}
\alias{GSE12470_eset}
\docType{data}
\title{Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis.}
\description{Gene expression profiles from 43 ovarian cancer tissues comprising 8 early stage and 35 advanced stage tissues were determined using oligonucleotide microarrays of 18,716 genes.  5 genes (ZEB2, CDH1, LTBP2, COL16A1 and ACTA2) were extracted as candidates for prognostic factors associated with progression-free survival.  The findings suggest that the expressions of epithelial-mesenchymal transition-related genes such as ZEB2 and CDH1 may play important roles in the invasion process of advanced-stage serous ovarian cancer.}
\usage{data(GSE12470_eset)}
\format{

Authors: Yoshihara K, Tajima A, Komata D, Yamamoto T, Kodama S, Fujiwara H, Suzuki M, Onishi Y, Hatae M, Sueyoshi K, Fujiwara H, Kudo Y, Inoue I, Tanaka K.

Lab: Department of Obstetrics and Gynecology, Niigata University Graduate School of Medical and Dental Sciences, Niigata, Japan.

Contact: tanaken@med.niigata-u.ac.jp

Title: Gene expression profiling of advanced-stage serous ovarian cancers distinguishes novel subclasses and implicates ZEB2 in tumor progression and prognosis.

Url: http://pubget.com/search?q=19486012

PubMedID: 19486012
}
\details{
Platform used: Agilent-012097 Human 1A Microarray (V2) G4110B (Feature Number version).

The data normalization was carried out using GeneSpring GX 7.3 (Agilent Technologies).  Processed data here are those provided by the original study authors.


assayData: 16788 features, 53 samples.

GEO_platform_accession: GSE12470.
platform_summary: GPL887.
biomart_ID: Agilent G4110B agilent_wholegenome.

---------------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------------

alt_sample_name:
Length = 53.
Class = character.
Mode = character.

sample_type: 
healthy = 10.  tumor = 43.

subtype: 
ser = 43. NA's = 10.

primarysite: 
ov = 53.

summarystage: 
early = 8. late = 35. NA's = 10.

tumorstage: 
2 = 8. 4 = 35. NA's = 10.

debulking: 
unknown = 53.
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE12470
}
\examples{
data(GSE12470_eset)
## maybe str(GSE12470_eset) ; plot(GSE12470_eset) ...
}
\keyword{datasets}
