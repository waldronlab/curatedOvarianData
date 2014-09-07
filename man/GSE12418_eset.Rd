\name{ GSE12418_eset }
\alias{ GSE12418_eset }
\docType{data}
\title{ Expression analysis of stage III serous ovarian adenocarcinoma distinguishes a sub-group of survivors. }
\description{ It is difficult to predict the clinical outcome for patients with ovarian cancer. However, the use of biomarkers as additional prognostic factors may improve the outcome for these patients. In order to find novel candidate biomarkers, differences in gene expressions were analysed in 54 stage III serous ovarian adenocarcinomas with oligonucleotide microarrays containing 27,000 unique probes. The microarray data was verified with quantitative real-time polymerase chain reaction for the genes TACC1, MUC5B and PRAME. Using hierarchical cluster analysis we detected a sub-group that included 60\% of the survivors. The gene expressions in tumours from patients in this sub-group of survivors were compared with the remaining tumours, and 204 genes were found to be differently expressed. We conclude that the sub-group of survivors might represent patients with favourable tumour biology and sensitivity to treatment. A selection of the 204 genes might be used as a predictive model to distinguish patients within and outside of this group. Alternative chemotherapy strategies could then be offered as first-line treatment, which may lead to improvements in the clinical outcome for these patients. }
\usage{data( GSE12418_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Partheen K, Levan K, Osterberg L, Horvath G.Expression analysis of stage III serous ovarian adenocarcinoma distinguishes a sub-group of survivors. Eur J Cancer. 2006 Nov; 42(16):2846-54. 
  Laboratory: Partheen, Horvath 2006 
  Contact information:  
  Title: Expression analysis of stage III serous ovarian adenocarcinoma distinguishes a sub-group of survivors. 
  URL:  
  PMIDs: 16996261 

  Abstract: A 177 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      SWEGENE H_v2.1.1_27k
   platform_shorttitle:     
      SWEGENE H_v2.1.1_27k
   platform_summary:     
      PartheenMetaData
   platform_manufacturer:     
      other
   platform_distribution:     
      non-commercial
   platform_accession:     
      GPL5886
   platform_technology:     
      spotted oligonucleotide

Preprocessing: default
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (12681 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 12681 features, 54 samples
Platform type: PartheenMetaData
Binary overall survival summary (definitions of long and short provided by study authors): 

 long short 
   20    34 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       54 character character 

sample_type: 
tumor 
   54 

histological_type: 
ser 
 54 

primarysite: 
ov 
54 

summarystage: 
late 
  54 

tumorstage: 
 3 
54 

substage: 
 b  c 
19 35 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  35.00   51.25   59.50   59.56   69.75   84.00 

pltx: 
 y 
54 

os_binary: 
 long short 
   20    34 

debulking: 
   optimal suboptimal 
        13         41 

uncurated_author_metadata: 
   Length     Class      Mode 
       54 character character 

}}
\keyword{datasets}
