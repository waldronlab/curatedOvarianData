\name{ GSE44104_eset }
\alias{ GSE44104_eset }
\docType{data}
\title{ COL11A1 promotes tumor progression and predicts poor clinical outcome in ovarian cancer. }
\description{ Biomarkers that predict disease progression might assist the development of better therapeutic strategies for aggressive cancers, such as ovarian cancer. Here, we investigated the role of collagen type XI alpha 1 (COL11A1) in cell invasiveness and tumor formation and the prognostic impact of COL11A1 expression in ovarian cancer. Microarray analysis suggested that COL11A1 is a disease progression-associated gene that is linked to ovarian cancer recurrence and poor survival. Small interference RNA-mediated specific reduction in COL11A1 protein levels suppressed the invasive ability and oncogenic potential of ovarian cancer cells and decreased tumor formation and lung colonization in mouse xenografts. A combination of experimental approaches, including real-time RT-PCR, casein zymography and chromatin immunoprecipitation (ChIP) assays, showed that COL11A1 knockdown attenuated MMP3 expression and suppressed binding of Ets-1 to its putative MMP3 promoter-binding site, suggesting that the Ets-1-MMP3 axis is upregulated by COL11A1. Transforming growth factor (TGF)-beta (TGF-??1) treatment triggers the activation of smad2 signaling cascades, leading to activation of COL11A1 and MMP3. Pharmacological inhibition of MMP3 abrogated the TGF-??1-triggered, COL11A1-dependent cell invasiveness. Furthermore, the NF-YA-binding site on the COL11A1 promoter was identified as the major determinant of TGF-??1-dependent COL11A1 activation. Analysis of 88 ovarian cancer patients indicated that high COL11A1 mRNA levels are associated with advanced disease stage. The 5-year recurrence-free and overall survival rates were significantly lower (P=0.006 and P=0.018, respectively) among patients with high expression levels of tissue COL11A1 mRNA compared with those with low expression. We conclude that COL11A1 may promote tumor aggressiveness via the TGF-??1-MMP3 axis and that COL11A1 expression can predict clinical outcome in ovarian cancer patients. }
\usage{data( GSE44104_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Wu Y, Chang T, Huang Y, Huang H, Chou C 
  Laboratory: Wu, Chou 2013 
  Contact information:  
  Title: COL11A1 promotes tumor progression and predicts poor clinical outcome in ovarian cancer. 
  URL:  
  PMIDs: 23934190 

  Abstract: A 260 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      Affymetrix HG-U133Plus2
   platform_summary:     
      hgu133plus2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL570
   platform_technology:     
      in situ oligonucleotide

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1BG-AS1 ... ZZZ3 (19816 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19816 features, 60 samples
Platform type: hgu133plus2
Binary overall survival summary (definitions of long and short provided by study authors): 

 long short 
   44    16 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       60 character character 

sample_type: 
tumor 
   60 

histological_type: 
clearcell      endo  mucinous       ser 
       12        11         9        28 

summarystage: 
early  late 
   25    35 

tumorstage: 
 1  2  3  4 
17  8 30  5 

recurrence_status: 
norecurrence   recurrence 
          40           20 

os_binary: 
 long short 
   44    16 

relapse_binary: 
 long short 
   40    20 

batch: 
2010-09-07 2010-09-08 2010-10-14 2010-12-10 2010-12-14 
        20          2         18         16          4 

uncurated_author_metadata: 
   Length     Class      Mode 
       60 character character 

}}
\keyword{datasets}
