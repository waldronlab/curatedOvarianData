\name{ GSE26193_eset }
\alias{ GSE26193_eset }
\docType{data}
\title{ miR-141 and miR-200a act on ovarian tumorigenesis by controlling oxidative stress response. }
\description{ Although there is evidence that redox regulation has an essential role in malignancies, its impact on tumor prognosis remains unclear. Here we show crosstalk between oxidative stress and the miR-200 family of microRNAs that affects tumorigenesis and chemosensitivity. miR-141 and miR-200a target p38?? and modulate the oxidative stress response. Enhanced expression of these microRNAs mimics p38?? deficiency and increases tumor growth in mouse models, but it also improves the response to chemotherapeutic agents. High-grade human ovarian adenocarcinomas that accumulate miR-200a have low concentrations of p38?? and an associated oxidative stress signature. The miR200a-dependent stress signature correlates with improved survival of patients in response to treatment. Therefore, the role of miR-200a in stress could be a predictive marker for clinical outcome in ovarian cancer. In addition, although oxidative stress promotes tumor growth, it also sensitizes tumors to treatment, which could account for the limited success of antioxidants in clinical trials. }
\usage{data( GSE26193_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Mateescu B, Batista L, Mariani O, Meyniel J, Cottu PH, Sastre-Garau X, Mechta-Grigoriou F 
  Laboratory: Mateescu, Mechta-Grigoriou 2011 
  Contact information:  
  Title: miR-141 and miR-200a act on ovarian tumorigenesis by controlling oxidative stress response. 
  URL:  
  PMIDs: 22101765 

  Abstract: A 149 word abstract is available. Use 'abstract' method.
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
assayData: 19816 features, 107 samples
Platform type: hgu133plus2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
 107.00  107.00  107.00   76.00    3.05    2.50    4.56 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      107 character character 

sample_type: 
tumor 
  107 

histological_type: 
clearcell      endo  mucinous     other       ser 
        6         8         8         6        79 

summarygrade: 
high  low 
  67   40 

summarystage: 
early  late 
   31    76 

tumorstage: 
 1  2  3  4 
20 11 59 17 

substage: 
   a    b    c NA's 
  16   12   62   17 

grade: 
 1  2  3 
 7 33 67 

days_to_tumor_recurrence: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    3.0   340.5   584.0  1108.0  1525.0  7386.0 

recurrence_status: 
norecurrence   recurrence 
          27           80 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      3     668    1096    1520    2220    7386 

vital_status: 
deceased   living 
      76       31 

batch: 
   Length     Class      Mode 
      107 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      107 character character 

}}
\keyword{datasets}
