\name{PMID15897565_eset}
\alias{PMID15897565_eset}
\docType{data}
\title{Patterns of gene expression that characterize long-term survival in advanced stage serous ovarian cancers.}
\description{A better understanding of the underlying biology of invasive serous ovarian cancer is critical for the development of early detection strategies and new therapeutics. The objective of this study was to define gene expression patterns associated with favorable survival.}
\usage{data(PMID15897565_eset)}
\format{

Authors: Berchuck A, Iversen ES, Lancaster JM, Pittman J, Luo J, Lee P, Murphy S, Dressman HK, Febbo PG, West M, Nevins JR, Marks JR.

Lab: Department of Obstetrics and Gynecology/Division of Gynecologic Oncology, Institute of Statistics and Decision Sciences, Center for Applied Genomics and Technology, Duke University Medical Center, Durham, North Carolina, USA.

Contact: berch001@mc.duke.edu

Title: Patterns of gene expression that characterize long-term survival in advanced stage serous ovarian cancers.

Url: http://pubget.com/search?q=15897565

PubMedID: 15897565

}
\details{
Platform used: HG-U133A.

frma normalization used from the frma bioconductor package.

assayData: 12858 features, 63 samples.

GEO_platform_accession: PMID15897565 
             
platform_summary: GPL96
             
biomart_ID: HG-U133A affy_hg_u133a

Binary overall survival summary (definitions of long and short provided by study authors):
(Long is when survival > 7 yrs); (Short is when survival < 3 years)

long = 24. short = 28. NA's = 11.

--------------------------------------------------------------------------------


Available sample meta-data: 

--------------------------------------------------------------------------------

alt_sample_name: 
Min. = 1761    1st Qu. = 1828   Median = 1907   Mean = 2001  3rd Qu. = 2032    Max. = 2536 

primarysite: 
ov = 63.

summarygrade: 
high = 25.  low = 37. NA's = 1.

summarystage: 
early = 11. late = 52.

tumorstage:
1 = 7. 2 = 4. 3 = 48. 4 = 4.

G: 
1 = 2. 2 = 35. 3 = 24. 4 = 1. NA's = 1.

age_at_initial_pathologic_diagnosis: 
Min. = 33.00  1st Qu. = 52.50  Median = 59.00  Mean = 59.21  3rd Qu. = 67.00   Max. = 79.00 

days_to_death: 
1092 = 28. 2555 = 24. NA's = 11.

os_binary: 
long = 24. short = 28. NA's = 11.

debulking: 
optimal = 24. suboptimal = 28. unknown = 11.

batch: 

2002-09-20 2002-10-23 2002-11-12 2002-12-16 2002-12-21 2003-01-03 2003-05-30 

        15          9         10          1          3         11         13 
        
2003-07-02 

         1 
}
\examples{
data(PMID15897565_eset)
## maybe str(PMID15897565_eset) ; plot(PMID15897565_eset) ...
#ISSUE: NO vital_status...

if(require(affy)){
  summary(PMID15897565_eset$os_binary)
}
}
\keyword{datasets}
