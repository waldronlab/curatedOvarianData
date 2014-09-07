\name{ GSE26712_eset }
\alias{ GSE26712_eset }
\docType{data}
\title{ A gene signature predicting for survival in suboptimally debulked patients with ovarian cancer. }
\description{ Despite the existence of morphologically indistinguishable disease, patients with advanced ovarian tumors display a broad range of survival end points. We hypothesize that gene expression profiling can identify a prognostic signature accounting for these distinct clinical outcomes. To resolve survival-associated loci, gene expression profiling was completed for an extensive set of 185 (90 optimal/95 suboptimal) primary ovarian tumors using the Affymetrix human U133A microarray. Cox regression analysis identified probe sets associated with survival in optimally and suboptimally debulked tumor sets at a P value of <0.01. Leave-one-out cross-validation was applied to each tumor cohort and confirmed by a permutation test. External validation was conducted by applying the gene signature to a publicly available array database of expression profiles of advanced stage suboptimally debulked tumors. The prognostic signature successfully classified the tumors according to survival for suboptimally (P = 0.0179) but not optimally debulked (P = 0.144) patients. The suboptimal gene signature was validated using the independent set of tumors (odds ratio, 8.75; P = 0.0146). To elucidate signaling events amenable to therapeutic intervention in suboptimally debulked patients, pathway analysis was completed for the top 57 survival-associated probe sets. For suboptimally debulked patients, confirmation of the predictive gene signature supports the existence of a clinically relevant predictor, as well as the possibility of novel therapeutic opportunities. Ultimately, the prognostic classifier defined for suboptimally debulked tumors may aid in the classification and enhancement of patient outcome for this high-risk population. }
\usage{data( GSE26712_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Bonome T, Levine DA, Shih J, Randonovich M, Pise-Masison CA, Bogomolniy F, Ozbun L, Brady J, Barrett JC, Boyd J, Birrer MJ: A Gene Signature Predicting for Survival in Suboptimally Debulked Patients with Ovarian Cancer. Cancer Research 2008, 68:5478 -5486. 
  Laboratory: Bonome, Birrer 2008 
  Contact information:  
  Title: A gene signature predicting for survival in suboptimally debulked patients with ovarian cancer. 
  URL:  
  PMIDs: 18593951 

  Abstract: A 238 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG-U133A] Affymetrix Human Genome U133A Array
   platform_shorttitle:     
      Affymetrix HG-U133A
   platform_summary:     
      hgu133a
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL96
   platform_technology:     
      in situ oligonucleotide

Preprocessing: frma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: A1CF A2M ... ZZZ3 (13104 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 13104 features, 195 samples
Platform type: hgu133a
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

   10 observations deleted due to missingness 
records   n.max n.start  events  median 0.95LCL 0.95UCL 
 185.00  185.00  185.00  129.00    3.83    3.24    4.83 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
      195 character character 

sample_type: 
healthy   tumor 
     10     185 

histological_type: 
 ser NA's 
 185   10 

primarysite: 
 ov 
195 

summarygrade: 
high NA's 
 185   10 

summarystage: 
late NA's 
 185   10 

tumorstage: 
   3    4 NA's 
 146   36   13 

substage: 
   b    c NA's 
   9  137   49 

age_at_initial_pathologic_diagnosis: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  26.00   52.00   63.00   61.54   70.00   84.00      13 

recurrence_status: 
norecurrence   recurrence 
          42          153 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   21.9   660.6  1164.0  1429.0  1880.0  4982.0      10 

vital_status: 
deceased   living     NA's 
     129       56       10 

debulking: 
   optimal suboptimal       NA's 
        90         95         10 

percent_normal_cells: 
20- 
195 

percent_stromal_cells: 
20- 
195 

percent_tumor_cells: 
80+ 
195 

batch: 
   Length     Class      Mode 
      195 character character 

uncurated_author_metadata: 
   Length     Class      Mode 
      195 character character 

}}
\keyword{datasets}
