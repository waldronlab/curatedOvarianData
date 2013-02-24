\name{ GSE18520_eset }
\alias{ GSE18520_eset }
\docType{data}
\title{ A gene signature predictive for outcome in advanced ovarian cancer identifies a survival factor: microfibril-associated glycoprotein 2. }
\description{ Advanced stage papillary serous tumors of the ovary are responsible for the majority of ovarian cancer deaths, yet the molecular determinants modulating patient survival are poorly characterized. Here, we identify and validate a prognostic gene expression signature correlating with survival in a series of microdissected serous ovarian tumors. Independent evaluation confirmed the association of a prognostic gene microfibril-associated glycoprotein 2 (MAGP2) with poor prognosis, whereas in vitro mechanistic analyses demonstrated its ability to prolong tumor cell survival and stimulate endothelial cell motility and survival via the alpha(V)beta(3) integrin receptor. Increased MAGP2 expression correlated with microvessel density suggesting a proangiogenic role in vivo. Thus, MAGP2 may serve as a survival-associated target. }
\usage{data( GSE18520_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Mok SC, Bonome T, Vathipadiekal V, Bell A, Johnson ME, Wong KK, Park DC, Hao K, Yip DK, Donninger H, Ozbun L, Samimi G, Brady J, Randonovich M, Pise-Masison CA, Barrett JC, Wong WH, Welch WR, Berkowitz RS, Birrer MJ.A gene signature predictive for outcome in advanced ovarian cancer identifies a survival factor: microfibril-associated glycoprotein 2. Cancer Cell. 2009 Dec 8; 16(6):521-32. 
  Laboratory: Mok, Birrer 2009 
  Contact information:  
  Title: A gene signature predictive for outcome in advanced ovarian cancer identifies a survival factor: microfibril-associated glycoprotein 2. 
  URL:  
  PMIDs: 19962670 

  Abstract: A 110 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
   platform_shorttitle:     
      Affymetrix HG-U133Plus2
   platform_summary:     
      hgu133plus2
   platform_manufacturer:     
      Affymetrix|Operon
   platform_distribution:     
      commercial|non-commercial
   platform_accession:     
      GPL570|GPL9216

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1BG A1CF ... ZZZ3 (19093 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 19093 features, 63 samples
Platform type: hgu133plus2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   10 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
  53.00   53.00   53.00   41.00    2.05    1.48    3.70 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  312.0   395.0   694.0   893.3  1040.0  2237.0 

sample_type: 
healthy   tumor 
     10      53 

histological_type: 
 ser NA's 
  53   10 

primarysite: 
ov 
63 

summarygrade: 
high NA's 
  53   10 

summarystage: 
late NA's 
  53   10 

tumorstage: 
   3 NA's 
  53   10 

grade: 
   3 NA's 
  53   10 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    150     450     630    1212    1440    4500      10 

vital_status: 
deceased   living     NA's 
      41       12       10 

debulking: 
optimal 
     63 

percent_normal_cells: 
 0 
63 

percent_stromal_cells: 
 0 
63 

percent_tumor_cells: 
100 
 63 

batch: 
2004-03-12 2004-04-08 2004-04-09 2004-07-20 2004-08-12 2004-08-13 2004-09-30 
        20          6          9         11         10          1          6 

uncurated_author_metadata: 
   Length     Class      Mode 
       63 character character 

}}
\keyword{datasets}
