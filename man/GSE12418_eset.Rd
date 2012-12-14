\name{GSE12418_eset}
\alias{GSE12418_eset}
\docType{data}
\title{Expression analysis of stage III serous ovarian adenocarcinoma distinguishes a sub-group of survivors.}
\description{In order to find novel candidate biomarkers to improve the outcome for patients with ovarian cancer, this study analyzed differences in gene expressions in 54 stage III serous ovarian adenocarcinomas with oligonucleotide microarrays containing 27,000 unique probes.  The microarray data was verified with quantitative real time polymerase chain reaction for the genes TACC1, MUC5B and PRAME.}
\usage{data(GSE12418_eset)}
\format{

Authors: Partheen K, Levan K, Osterberg L, Horvath G.

Lab: Department of Oncology, Goteborg University, SE-413 45 Goteborg, Sweden.

Contact: karolina.partheen@oncology.gu.se

Title: Expression analysis of stage III serous ovarian adenocarcinoma distinguishes a sub-group of survivors.

Url: http://pubget.com/search?q=16996261

PubMedID: 16996261.
}
\details{
Platform used: SWEGENE H_v2.1.1_27k.

loess normalization was used.  Processed data here are those provided by the original study authors.

assayData: 12633 features, 54 samples.

GEO_platform_accession : GSE12418.

platform_summary: GPL5886.

biomart_ID: Swegene

Binary overall survival summary (definitions of long and short provided by study authors): 

long = 20. short = 34.

---------------------------------------------------------------------------------


Available sample meta-data: 

--------------------------------------------------------------------------------- 

alt_sample_name: 
Length = 54.
Class = character.
Mode = character.

sample_type: 
tumor = 54.

subtype: 
ser = 54.

primarysite: 
ov = 54. 

summarystage: 
late = 54.

tumorstage: 
3 = 54. 

substage: 
b = 19. c = 35. 

age_at_initial_pathologic_diagnosis: 
Min. = 35.00   1st Qu. = 51.25  Median = 59.50   Mean = 59.56  3rd Qu. = 69.75   Max = 84.00 

pltx: 
y = 54. 

os_binary: 
long = 20. short = 34. 
   
debulking: 
optimal = 13. suboptimal = 41.
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12418
}
\examples{
data(GSE12418_eset)
## maybe str(GSE12418_eset) ; plot(GSE12418_eset) ...
if(require(affy)){
  summary(GSE12418_eset$os_binary)
}
}
\keyword{datasets}
