\name{E.MTAB.386_eset}
\alias{E.MTAB.386_eset}
\docType{data}
\title{Angiogenic mRNA and microRNA gene expression signatures predict
  overall survival in serous ovarian cancer (mRNA data only).}
\description{Ovarian cancer is the fifth leading cause of cancer death for women in the U.S. and the seventh most fatal worldwide. Although ovarian cancer is notable for its initial sensitivity to platinum-based therapies, the vast majority of women eventually recurs and succumbs to increasingly platinum-resistant disease. Modern, targeted cancer drugs intervene in cell signaling pathways by compensating for their deregulation, and identifying these key mechanisms and pathways would greatly advance our ability to treat disease. In order to shed light on the molecular diversity of ovarian cancer, we performed comprehensive transcriptional profiling on 129 high grade, late stage serous ovarian cancers. We implemented a novel, re-sampling based version of the ISIS class discovery algorithm (rISIS: robust ISIS) and applied it to the entire set of ovarian cancer transcriptional profiles. rISIS identified a novel stratification of this disease into two groups with significantly different overall survival. Gene set enrichment analysis found strong support for the stratification by extracellular matrix, cell adhesion, and angiogenesis genes. Application of this 'angiogenesis' signature to independent, published ovarian cancer gene expression data confirms its prognostic potential. Additional support for this stratification is provided by micro-RNA expression profiles which exhibit statistically significant expression differences between the groups, and additional mechanistic analyses have allowed development of hypotheses relevant to directed therapeutic intervention for specific subclasses of the disease. In particular, the subgroup stratification we discovered may be relevant for identifying which patients may be best suited for anti-angiogenic therapies that are now being tested in clinical trials.
}
\usage{data(E.MTAB.386_eset)}
\format{

Authors: Stefan Bentink, John Quackenbush, Jian-Bing Fan, Michelle S. Hirsch, Kristina Holton, Renee Rubio, Craig April, Jing Chen, Eliza Wickham-Garcia, Joyce Liu, Aedin Culhane, Ronny Drapkin, and Ursula A. Matulonis.

Contact: John Quackenbush, Stefan Bentink <bentink@jimmy.harvard.edu>, Ursula A Matulonis <ursula_matulonis@dfci.harvard.edu>

Title: Angiogenic mRNA and microRNA gene expression signatures predict overall survival in serous ovarian cancer.

Url: http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-386

PubMedID: 22348002
}
\details{
Platforms used: A-MEXP-1811 - Illumina Human microRNA MI_V1_R0_XS0000122-MAP & A-MEXP-931 - Illumina HumanRef-8 v2 Expression BeadChip.


Overall survival time-to-event summary (in years): 

Call: survfit(formula = Surv(time, cens) ~ -1)

---------------------------------------------------------------------------------


Available sample meta-data: 

---------------------------------------------------------------------------------

alt_sample_name: 
Length = 141.
Class = character.
Mode = character.

sample_type:
tumor = 141.

subtype:
ser = 141.
   
primarysite: 
ov = 141.

summarystage: 
late = 140. NA's = 1.

tumorstage:
3 = 118. 4 = 22. NA's = 1.

substage: 
a = 4.    b = 14. c = 101. NA's = 22.

age_at_initial_pathologic_diagnosis: 
Min. = 21.21  Median = 66.42  Mean = 60.84  Max. = 95.13  
 
days_to_death: 
Min. = 3.9  Median =  920.1  Mean = 1039.36  Max. = 2724 

vital_status: 
deceased = 81. living = 60. 

debulking: 
optimal = 34. suboptimal = 107. 
}
\source{
http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-386
}
\examples{
data(E.MTAB.386_eset)
## maybe str(GSE12418_eset) ; plot(GSE12418_eset) ...
if(require(affy)){
  summary(E.MTAB.386_eset$os_binary)
}
}
\keyword{datasets}
