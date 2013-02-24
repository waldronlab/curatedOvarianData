\name{ GSE19829.GPL8300_eset }
\alias{ GSE19829.GPL8300_eset }
\docType{data}
\title{ Gene expression profile of BRCAness that correlates with responsiveness to chemotherapy and with outcome in patients with epithelial ovarian cancer. }
\description{ To define a gene expression profile of BRCAness that correlates with chemotherapy response and outcome in epithelial ovarian cancer (EOC).A publicly available microarray data set including 61 patients with EOC with either sporadic disease or BRCA(1/2) germline mutations was used for development of the BRCAness profile. Correlation with platinum responsiveness was assessed in platinum-sensitive and platinum-resistant tumor biopsy specimens from six patients with BRCA germline mutations. Association with poly-ADP ribose polymerase (PARP) inhibitor responsiveness and with radiation-induced RAD51 foci formation (a surrogate of homologous recombination) was assessed in Capan-1 cell line clones. The BRCAness profile was validated in 70 patients enriched for sporadic disease to assess its association with outcome.The BRCAness profile accurately predicted platinum responsiveness in eight out of 10 patient-derived tumor specimens, and between PARP-inhibitor sensitivity and resistance in four out of four Capan-1 clones. [corrected] When applied to the 70 patients with sporadic disease, patients with the BRCA-like (BL) profile had improved disease-free survival (34 months v 15 months; log-rank P = .013) and overall survival (72 months v 41 months; log-rank P = .006) compared with patients with a non-BRCA-like (NBL) profile, respectively. The BRCAness profile maintained independent prognostic value in multivariate analysis, which controlled for other known clinical prognostic factors.The BRCAness profile correlates with responsiveness to platinum and PARP inhibitors and identifies a subset of sporadic patients with improved outcome. Additional evaluation of this profile as a predictive tool in patients with sporadic EOC is warranted. }
\usage{data( GSE19829.GPL8300_eset )}
\format{
\preformatted{
experimentData(eset):
Experiment data
  Experimenter name: Konstantinopoulos PA, Spentzos D, Karlan BY, Taniguchi T et al. Gene expression profile of BRCAness that correlates with responsiveness to chemotherapy and with outcome in patients with epithelial ovarian cancer. J Clin Oncol 2010 Aug 1;28(22):3555-61. 
  Laboratory: Konstantinopoulos, Cannistra 2010 hgu95 
  Contact information:  
  Title: Gene expression profile of BRCAness that correlates with responsiveness to chemotherapy and with outcome in patients with epithelial ovarian cancer. 
  URL:  
  PMIDs: 20547991 

  Abstract: A 241 word abstract is available. Use 'abstract' method.
  Information is available on: preprocessing 
  notes:
   platform_title:     
      [HG_U95Av2] Affymetrix Human Genome U95 Version 2 Array
   platform_shorttitle:     
      Affymetrix HG_U95Av2
   platform_summary:     
      hgu95av2
   platform_manufacturer:     
      Affymetrix
   platform_distribution:     
      commercial
   platform_accession:     
      GPL570|GPL8300

Preprocessing: rma
featureData(eset):
An object of class 'AnnotatedDataFrame'
  featureNames: AADAC AAK1 ... ZZZ3 (8957 total)
  varLabels: probeset gene
  varMetadata: labelDescription

}}
\details{
\preformatted{
assayData: 8957 features, 42 samples
Platform type: hgu95av2
Overall survival time-to-event summary (in years):
Call: survfit(formula = Surv(time, cens) ~ -1)

records   n.max n.start  events  median 0.95LCL 0.95UCL 
  42.00   42.00   42.00   23.00    3.78    2.79      NA 

--------------------------- 
Available sample meta-data: 
--------------------------- 

alt_sample_name: 
   Length     Class      Mode 
       42 character character 

sample_type: 
tumor 
   42 

primarysite: 
ov 
42 

days_to_death: 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   30.0   727.5  1155.0  1089.0  1485.0  2040.0 

vital_status: 
deceased   living 
      23       19 

batch: 
2001-09-14 2001-12-14 2002-08-20 2003-09-09 2003-09-18 
         7          4         14         13          4 

uncurated_author_metadata: 
   Length     Class      Mode 
       42 character character 

}}
\keyword{datasets}
