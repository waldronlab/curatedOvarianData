\name{GSE2109_eset}
\alias{GSE2109_eset}
\docType{data}
\title{Expression Project for Oncology (expO)}
\description{This project seeks to integrate longitudinal clinical annotation with gene expression data for a unique and powerful portrait of human malignancies, providing critical perspective on diagnositic markers, prognostic indicators and therapeutic targets.  Tissue samples were procured under standard conditions and gene expression analyses were performed on a clinically annotated set of tumor samples.
}
\usage{data(GSE2109_eset)}
\format{

Title: The Expression Project for Oncology (EXPO)

Contact: ecurley@intgen.org

Url: http://www.intgen.org/expo/

}
\details{
Platform used: [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array.

frma normalization used from the frma bioconductor package.

assayData: 18708 features, 204 samples

GEO_platform_accession: GSE2109

platform_summary: GPL570          

biomart_ID: HG-U133_Plus_2 affy_hg_u133_plus_2

------------------------------------------------------------------------------- 


Available sample meta-data: 

-------------------------------------------------------------------------------


alt_sample_name: 
Length= 204.
Class= character.
Mode= character.

sample_type: 
tumor = 204. 
  
subtype: 
Length = 204.
Class = character.
Mode = character.
  
primarysite: 
other= 23.  ov= 178.  unknown= 3. 
     
summarygrade: 
high= 91.  low= 31.  NA's= 82. 
  
summarystage: 
early= 37.   late= 85.   NA's= 82. 

T:
1= 21.    2= 16.    3= 85.    NA's= 82. 

substage: 
a = 17.    b= 22.    c= 79.    NA's = 86. 
  
G: 
1 = 11.    2= 20.    3= 83.    4= 8.   NA's= 82. 
  
N: 
0 = 60.    1 = 43.     2= 2.   NA's= 99.
  
M: 
0 = 92.    1 = 18.    NA's = 94.

age_at_initial_pathologic_diagnosis: 
Min.= 25.00   1st Qu.= 45.00  Median= 55.00    Mean= 58.82    3rd Qu.= 65.00    Max.= 85.00
  
debulking: 
unknown = 204 
}
\source{
http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE2109
}
\examples{
data(GSE2109_eset)
## maybe str(GSE2109_eset) ; plot(GSE2109_eset) ...
}
\keyword{datasets}
