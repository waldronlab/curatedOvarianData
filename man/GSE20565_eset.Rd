\name{ GSE20565_eset }
\alias{ GSE20565_eset }
\docType{data}
\title{ A genomic and transcriptomic approach for a differential diagnosis between primary and secondary ovarian carcinomas in patients with a previous history of breast cancer. }
\description{ The distinction between primary and secondary ovarian tumors may be challenging for pathologists. The purpose of the present work was to develop genomic and transcriptomic tools to further refine the pathological diagnosis of ovarian tumors after a previous history of breast cancer.Sixteen paired breast-ovary tumors from patients with a former diagnosis of breast cancer were collected. The genomic profiles of paired tumors were analyzed using the Affymetrix GeneChip Mapping 50 K Xba Array or Genome-Wide Human SNP Array 6.0 (for one pair), and the data were normalized with ITALICS (ITerative and Alternative normaLIzation and Copy number calling for affymetrix Snp arrays) algorithm or Partek Genomic Suite, respectively. The transcriptome of paired samples was analyzed using Affymetrix GeneChip Human Genome U133 Plus 2.0 Arrays, and the data were normalized with gc-Robust Multi-array Average (gcRMA) algorithm. A hierarchical clustering of these samples was performed, combined with a dataset of well-identified primary and secondary ovarian tumors.In 12 of the 16 paired tumors analyzed, the comparison of genomic profiles confirmed the pathological diagnosis of primary ovarian tumor (n = 5) or metastasis of breast cancer (n = 7). Among four cases with uncertain pathological diagnosis, genomic profiles were clearly distinct between the ovarian and breast tumors in two pairs, thus indicating primary ovarian carcinomas, and showed common patterns in the two others, indicating metastases from breast cancer. In all pairs, the result of the transcriptomic analysis was concordant with that of the genomic analysis.In patients with ovarian carcinoma and a previous history of breast cancer, SNP array analysis can be used to distinguish primary and secondary ovarian tumors. Transcriptomic analysis may be used when primary breast tissue specimen is not available. }
\usage{data( GSE20565_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Meyniel JP, Cottu PH, Decraene C, Stern MH, Couturier J, Lebigot I, Nicolas A, Weber N, Fourchotte V, Alran S, Rapinat A, Gentien D, Roman-Roman S, Mignot L, Sastre-Garau X.A genomic and transcriptomic approach for a differential diagnosis between primary and secondary ovarian carcinomas in patients with a previous history of breast cancer. BMC Cancer. 2010 May 21; 10:222. 
  Laboratory: Meyniel, Sastre-Garau 2010 
  Contact information:  
  Title: A genomic and transcriptomic approach for a differential diagnosis between primary and secondary ovarian carcinomas in patients with a previous history of breast cancer. 
  URL:  
  PMIDs: 20492709 

  Abstract: A 277 word abstract is available. Use 'abstract' method.
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
assayData: 19816 features, 140 samples
Platform type: hgu133plus2
--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      140 character character 

sample_type: 
tumor 
  140 

histological_type: 
   Length     Class      Mode 
      140 character character 

primarysite: 
other    ov 
   44    96 

summarygrade: 
high  low NA's 
  63   33   44 

summarystage: 
early  late  NA's 
   27    67    46 

tumorstage: 
   1    2    3    4 NA's 
  18    9   52   15   46 

substage: 
   a    b    c NA's 
  14   10   55   61 

grade: 
   1    2    3 NA's 
   6   27   63   44 

batch: 
   Length     Class      Mode 
      140 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      140 character character 

}}
\keyword{datasets}
